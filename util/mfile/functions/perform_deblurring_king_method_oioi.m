function [im_king, im] = perform_deblurring_king_method_oioi(im, kx, ky, gx, gy, gz, X_, Y_, Z_, field_of_view_mm, DCS_offset, gamma, B0, dt, A)
% Written by Nam Gyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 12/27/2020, Last modified: 12/27/2020
% Modified by Ye Tian, 06/07/2021

%% Get data dimension
X = X_;
Y = Y_;
Z = Z_;

[Nk, Ni] = size(kx);
N1      = size(im, 1);
N2      = size(im, 2);
Nframe  = size(im, 3);
Necho   = size(im, 4);

im = fliplr(rot90(im, -1));

%% Index for spiral out in out in
idx1 = 1:round(Nk/4);
idx2 = round(Nk/4)+1:Nk/2;
idx3 = Nk/2+1:Nk/2+round(Nk/4);
idx4 = Nk/2+round(Nk/4)+1:Nk;

%% Calculate the transformation matrix from RCS to DCS
%--------------------------------------------------------------------------
% From the paper: xyz in the physical coordinates
% [x]   [a1 a2 a3][X]
% [y] = [a4 a5 a6][Y]
% [z]   [a7 a8 a9][Z]
%
% [gx]   [a1 a2 a3][GX]
% [gy] = [a4 a5 a6][GY]
% [gz]   [a7 a8 a9][GZ]
%
% For Siemens:
% [x]                                       [r]
% [y] = R_pcs2dcs * R_gcs2pcs * R_rcs2gcs * [c] + DCS_offset
% [z]                                       [s]
%
% xyz = A * rcs + A * A.' * DCS_offset
%     = A * (rcs + A.' * DCS_offset)
%     = A * XYZ
%
% [gx]                                       [gr]
% [gy] = R_pcs2dcs * R_gcs2pcs * R_rcs2gcs * [gc]
% [gz]                                       [gs]
%--------------------------------------------------------------------------
% A = rotMatrixPCSToDCS * rotMatrixGCSToPCS * rotMatrixRCSToGCS;

a1 = A(1,1); a2 = A(1,2); a3 = A(1,3);
a4 = A(2,1); a5 = A(2,2); a6 = A(2,3);
a7 = A(3,1); a8 = A(3,2); a9 = A(3,3);

%% Calculate constants (F1 to F6)
F1 =  1 / 4 * (a1^2 + a4^2) * (a7^2 + a8^2) + a7^2 * (a2^2 + a5^2) - a7 * a8 * (a1 * a2 + a4 * a5);
F2 =  1 / 4 * (a2^2 + a5^2) * (a7^2 + a8^2) + a8^2 * (a1^2 + a4^2) - a7 * a8 * (a1 * a2 + a4 * a5);
F3 =  1 / 4 * (a3^2 + a6^2) * (a7^2 + a8^2) + a9^2 * (a1^2 + a2^2 + a4^2 + a5^2) - a7 * a9 * (a1 * a3 + a4 * a6) - a8 * a9 * (a2 * a3 + a5 * a6);
F4 =  1 / 2 * (a2 * a3 + a5 * a6) * (a7^2 - a8^2) + a8 * a9 * (2 * a1^2 + a2^2 + 2 * a4^2 + a5^2) - a7 * a8 * (a1 * a3 + a4 * a6) - a7 * a9 * (a1 * a2 + a4 * a5);
F5 =  1 / 2 * (a1 * a3 + a4 * a6) * (a8^2 - a7^2) + a7 * a9 * (a1^2 + 2 * a2^2 + a4^2 + 2 * a5^2) - a7 * a8 * (a2 * a3 + a5 * a6) - a8 * a9 * (a1 * a2 + a4 * a5);
F6 = -1 / 2 * (a1 * a2 + a4 * a5) * (a7^2 + a8^2) + a7 * a8 * (a1^2 + a2^2 + a4^2 + a5^2);

%% Perform through-plane correction
% start_time = tic; fprintf('Performing through-plane correction... ');
%--------------------------------------------------------------------------
% Calculate a scaled concomitant field time parameter tc(t) [sec]
%--------------------------------------------------------------------------
% g0 = sqrt(GX.^2 + GY.^2); % Nk x Ni [G/cm]
g0 = sqrt(gx.^2 + gy.^2);
gm = max(g0(:)); % [G/cm]
tc = 1 / gm^2 * cumsum(g0.^2 * dt); % 1 / [G/cm]^2 * [G/cm]^2 * [sec] => [sec]

%--------------------------------------------------------------------------
% Calculate the time-independent frequency offset fc
%--------------------------------------------------------------------------
% XYZ_offset = A * DCS_offset;
% [rad/sec/T] / [2pi rad/cycle] * [m]^2 / [T] * ([G/cm] * [T/1e4G] * [1e2cm/m])^2 => [cycle/sec]
% fc = gamma / (2 * pi) * (gm * 1e-2)^2 / (4 * B0) * F3 * XYZ_offset(3)^2; % [Hz]

%--------------------------------------------------------------------------
% Calculate the through-plane concomitant field phase
%--------------------------------------------------------------------------
% phi_c = 2 * pi * fc * tc; % [rad/cycle] * [Hz] * [sec] => [rad]

%--------------------------------------------------------------------------
% Demodulation before gridding
%--------------------------------------------------------------------------
% kspace_demod = bsxfun(@times, exp(1j * phi_c), kspace);
% kspace_demod = kspace;
% fprintf('done! (%6.4f sec)\n', toc(start_time));

%% Perform NUFFT reconstruction
%----------------------------------------------------------------------
% Apply the adjoint of 2D NUFFT
%----------------------------------------------------------------------
% imc_nufft = NUFFT.NUFFT_adj(kspace_demod, N);
% imc_nufft = fliplr(rot90(imc_nufft, -1));

%% IFFT to k-space (k-space <=> image-space)
kspace_cart = fftshift2(fft2(fftshift2(im)));

%% Perform in-plane correction (frequency-segmented deblurring)
L = 110; % number of frequency bins

%--------------------------------------------------------------------------
% Calculate the frequency offsets for concomitant fields
%--------------------------------------------------------------------------
fc_XYZ = gamma / (2 * pi) * (gm * 1e-2)^2 / (4 * B0) * (F1 * X.^2 + F2 * Y.^2 + F3 * Z.^2 + F4 * Y .* Z + F5 * X .* Z + F6 * X .* Y); % [Hz]

%--------------------------------------------------------------------------
% Estimate the range of frequency offsets
% [rad/sec/T] / [2pi rad/cycle] * ([G/cm] * [T/1e4G] * [1e2cm/m])^2 / [T] * [m]^2 => [cycle/sec]
%--------------------------------------------------------------------------
D = max(field_of_view_mm) * 1e-3; % [mm] * [m/1e3mm] => [m]
index = (X.^2 + Y.^2 < (D / 2).^2);
df_max = max(fc_XYZ(index)); % [Hz]
L = round(df_max);
interval = df_max / (L - 1); % [Hz]
df_range = (0:L-1).' * interval; % [Hz]
df_range = unique(df_range);
L = length(df_range);

%--------------------------------------------------------------------------
% Perform 2D interpolation on tc(t) [sec]
%--------------------------------------------------------------------------
KX_range = (-floor(N1/2):ceil(N1/2)-1).' * 1 / N1; % [-0.5,0.5]
KY_range = (-floor(N2/2):ceil(N2/2)-1).' * 1 / N2; % [-0.5,0.5]
[KY_grid,KX_grid,~] = ndgrid(KX_range, KY_range, 0);

F1 = scatteredInterpolant(double(vec(kx(idx1,:))), double(vec(ky(idx1,:))), double(vec(tc(idx1,:))), 'linear', 'linear');
tc_grid_1 = F1(KX_grid, KY_grid);

F2 = scatteredInterpolant(double(vec(kx(idx2,:))), double(vec(ky(idx2,:))), double(vec(tc(idx2,:))), 'linear', 'linear');
tc_grid_2 = F2(KX_grid, KY_grid);

F3 = scatteredInterpolant(double(vec(kx(idx3,:))), double(vec(ky(idx3,:))), double(vec(tc(idx3,:))), 'linear', 'linear');
tc_grid_3 = F3(KX_grid, KY_grid);

F4 = scatteredInterpolant(double(vec(kx(idx4,:))), double(vec(ky(idx4,:))), double(vec(tc(idx4,:))), 'linear', 'linear');
tc_grid_4 = F4(KX_grid, KY_grid);
%--------------------------------------------------------------------------
% Perform frequency-segmented deblurring
%--------------------------------------------------------------------------
index = fc_XYZ - permute(df_range, [2, 3, 1]);
index = abs(index);
[~, index] = min(index, [], 3);
im_king = zeros([N1, N2, Nframe, Necho], 'single');
im_fs = zeros([N1, N2, Nframe, Necho], 'single');
frequency_bins = unique(index(:))';

tc_grid_1 = gpuArray(tc_grid_1);
tc_grid_2 = gpuArray(tc_grid_2);
tc_grid_3 = gpuArray(tc_grid_3);
tc_grid_4 = gpuArray(tc_grid_4);
df_range = gpuArray(df_range);
kspace_cart = gpuArray(kspace_cart);
im_king = gpuArray(im_king);
im_fs = gpuArray(im_fs);

tic
for idx1 = frequency_bins
    %----------------------------------------------------------------------
    % Demodulation after gridding
    %----------------------------------------------------------------------
    phi_c = 2 * pi * df_range(idx1) * tc_grid_1; % [rad/cycle] * [Hz] * [sec] => [rad]
    kspace_fs_gridded = bsxfun(@times, exp(1j * phi_c), kspace_cart(:, :, :, 1));

    %----------------------------------------------------------------------
    % iFFT to image-space (k-space <=> image-space)
    %----------------------------------------------------------------------
    im_fs(:, :, :, 1) = fftshift2(ifft2(fftshift2(kspace_fs_gridded)));
    
    phi_c = 2 * pi * df_range(idx1) * tc_grid_2; % [rad/cycle] * [Hz] * [sec] => [rad]
    kspace_fs_gridded = bsxfun(@times, exp(1j * phi_c), kspace_cart(:, :, :, 2));
    im_fs(:, :, :, 2) = fftshift2(ifft2(fftshift2(kspace_fs_gridded)));
    
    phi_c = 2 * pi * df_range(idx1) * tc_grid_3; % [rad/cycle] * [Hz] * [sec] => [rad]
    kspace_fs_gridded = bsxfun(@times, exp(1j * phi_c), kspace_cart(:, :, :, 3));
    im_fs(:, :, :, 3) = fftshift2(ifft2(fftshift2(kspace_fs_gridded)));
    
    phi_c = 2 * pi * df_range(idx1) * tc_grid_4; % [rad/cycle] * [Hz] * [sec] => [rad]
    kspace_fs_gridded = bsxfun(@times, exp(1j * phi_c), kspace_cart(:, :, :, 4));
    im_fs(:, :, :, 4) = fftshift2(ifft2(fftshift2(kspace_fs_gridded)));
    
    im_king = im_king + im_fs .* (index == idx1);
    
    if mod(find(idx1 == frequency_bins), 100) == 0
        fprintf('progress: %02.f%%, ',   find(idx1 == frequency_bins) / length(frequency_bins) * 100)
        toc
        tic
    end
end

im_king = gather(im_king);
% if L == 1
%     im_king = im_fs;
%     imageMRI([sos(im), sos(im_king), sos(im) - sos(im_king)]);
%     return;
% end

%--------------------------------------------------------------------------
% Perform pixel-dependent interpolation (nearest neighbor)
%--------------------------------------------------------------------------
% index = fc_XYZ - permute(df_range, [2, 3, 1]);
% index = abs(index);
% [~, index] = min(index, [], 3);
% im_king = zeros([N1, N2, Nframe, Necho], 'single');
% for idxc = 1:Nframe
%     for idx2 = 1:N2
%         for idx1 = 1:N1
%             im_king(idx1, idx2, idxc, :) = im_fs(idx1, idx2, index(idx1, idx2), idxc, :);
%         end
%     end
% end

% imageMRI([im(:,:,round(Nframe/2),:); im_king(:,:,round(Nframe/2),:); abs(im(:,:,round(Nframe/2),:)) - abs(im_king(:,:,round(Nframe/2),:))]);

end