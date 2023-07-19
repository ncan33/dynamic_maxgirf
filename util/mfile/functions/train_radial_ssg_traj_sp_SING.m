function kSpace_cart = train_radial_ssg_traj_sp_SING(im_sb, theta, kxs, kys, kxt, kyt, phase_mod_idx, kSpace_mb)

fprintf([repmat('-', [1, 75]), '\n'])
fprintf(sprintf('begin trajectory specified radial slice GRAPPA calibration: \n'));
fprintf([repmat('-', [1, 75]), '\n'])

kernel_n_max = 25;
max_d = 5;
% lambda = 0;
lambda = 0.001;
kb_size = [2, 2];

% size
[sx, ~, nof_calib, nc, nsms] = size(im_sb);
kxs = kxs(:);
kys = kys(:);
kxt = kxt(:);
kyt = kyt(:);
n_target = length(kxt);

% single band cartesian k-space
N = NUFFT.init(kxs, kys, 1, kb_size, sx, sx);
kSpace_cart_calib = fft2(im_sb.*N.Apodizer); clear im_sb

% shift kernel
shift_kx = -10:2:10;
shift_ky = -10:2:10;
[shift_kx, shift_ky] = meshgrid(shift_kx, shift_ky);
shift_kx = shift_kx(:)';
shift_ky = shift_ky(:)';
nshift = length(shift_kx);

% phase modulation
phase_mod_idx = repmat(phase_mod_idx, [sx, 1]);
phase_mod_idx = phase_mod_idx(:);
phase_mod_ = exp(-1i*2*pi*phase_mod_idx/3);
phase_mod_ = cat(2, ones(size(phase_mod_)), phase_mod_, conj(phase_mod_));

% init
kSpace_cart = zeros([sx, sx, nc, nsms], 'like', kSpace_cart_calib);

% display some text
text_display = sprintf('number of calibration time frames: %4.0f\nnumber of coils: %22.0f\nnumber of slices: %21.0f\n\n', nof_calib, nc, nsms);
fprintf(text_display)

kernel = cell(n_target, 1);

kx_all_ = cell(n_target, 1);
ky_all_ = cell(n_target, 1);
nsource_ = zeros(n_target, 1);
phase_mod_source_ = cell(n_target, 1);
source_idx_all = cell(n_target, 1);
for ipoint = 1:n_target
    % find n points that are mostly close
    kx_target = kxt(ipoint);
    ky_target = kyt(ipoint);
    
    d = sqrt((kx_target - kxs).^2 + (ky_target - kys).^2);
    source_idx = find(d <= max_d);
    nsource = length(source_idx);
    if nsource > kernel_n_max
        random_pick = randperm(nsource);
        random_pick = random_pick(1:kernel_n_max);
        source_idx = source_idx(random_pick);
        nsource = kernel_n_max;
    elseif nsource == 0
        source_idx = find(min(d));
        nsource = length(source_idx);
    end
    
    phase_mod_source = phase_mod_(source_idx, :);
    phase_mod_source = reshape(phase_mod_source, [nsource, 1, 1, 1, nsms]);
    
    kx_all = [kx_target; kxs(source_idx)];
    ky_all = [ky_target; kys(source_idx)];
    
    % put the kernel to bottow left cornor in calibration area
    % and add all the shifts
    kx_all_{ipoint} = vec(kx_all + shift_kx);
    ky_all_{ipoint} = vec(ky_all + shift_ky);
    nsource_(ipoint) = nsource;
    phase_mod_source_{ipoint} = phase_mod_source;
    source_idx_all{ipoint} = source_idx;
    
    % test
%     n_point_all = (nsource + 1) * nshift;
%     kspace_calib = zeros([n_point_all, nc, nof_calib, nsms], 'like', kSpace_cart_calib);
%     
%     N = NUFFT.init(vec(kx_all + shift_kx), vec(ky_all + shift_ky), 1, kb_size, sx ,sx);
%     for itime = 1:nof_calib
%         kspace_calib(:, :, itime, :) = NUFFT.cart2rad(kSpace_cart_calib(:, :, itime, :, :), N);
%     end
%     kspace_calib = reshape(kspace_calib, [nsource + 1, nshift, nc, nof_calib, nsms]);
%     source = kspace_calib(2:end, :, :, :, :);
%     source = sum(source .* phase_mod_source, 5);
%     source = permute(source, [2, 4, 1, 3]);
%     source = reshape(source, [nshift*nof_calib, nsource * nc]);
%     scale = max(abs(source(:)));
%     source = source/scale;
%     
%     target = kspace_calib(1, :, :, :, :);
%     target = permute(target, [2, 4, 3, 5, 1]);
%     target = reshape(target, [nshift*nof_calib, nc*nsms])/scale;
%     
%     % Tikhonov regularizaiton
%     kernel{ipoint} = ((source' * source + lambda * eye(nsource*nc))^-1 * source' * target).';
end

npoint_per_seg = 300;
nseg = ceil(n_target/npoint_per_seg);

ipoint_total = 0;
for iseg = 1:nseg
    tic
    fprintf(sprintf('%g / %g ', iseg, nseg))
    if iseg < nseg
        idx_seg = (1:npoint_per_seg) + (iseg - 1)*npoint_per_seg;
    else
        idx_seg = npoint_per_seg * (nseg - 1) + 1 : n_target;
    end
    
    kx_all = cat(1, kx_all_{idx_seg});
    ky_all = cat(1, ky_all_{idx_seg});
    % interpolate calibration data
    n_point_all = length(kx_all);
    kspace_calib = zeros([n_point_all, nc, nof_calib, nsms], 'like', kSpace_cart_calib);
    
    N = NUFFT.init(kx_all, ky_all, 1, kb_size, sx ,sx);
    for itime = 1:nof_calib
        kspace_calib(:, :, itime, :) = NUFFT.cart2rad(kSpace_cart_calib(:, :, itime, :, :), N);
    end
    
    for ipoint = 1:length(idx_seg)
        % estimate GRAPPA kernel
        idx_point = idx_seg(ipoint);
%         if nsource_(idx_point)>0
            idx_points = (1:(nsource_(idx_point)+1)*nshift) + (sum(nsource_(1:idx_point-1)) + (idx_point - 1))*nshift - ipoint_total;
            kspace_calib_temp = reshape(kspace_calib(idx_points, :, :, :), [nsource_(idx_point) + 1, nshift, nc, nof_calib, nsms]);
            source = kspace_calib_temp(2:end, :, :, :, :);
            source = sum(source .* phase_mod_source_{idx_point}, 5);
            source = permute(source, [2, 4, 1, 3]);
            source = reshape(source, [nshift*nof_calib, nsource_(idx_point) * nc]);
            scale = max(abs(source(:)));
            source = source/scale;
            
            target = kspace_calib_temp(1, :, :, :, :);
            target = permute(target, [2, 4, 3, 5, 1]);
            target = reshape(target, [nshift*nof_calib, nc*nsms])/scale;
            
            % Tikhonov regularizaiton
            kernel{idx_point} = ((source' * source + lambda * eye(nsource_(idx_point)*nc))^-1 * source' * target).';
%         end
    end
    ipoint_total = ipoint_total + idx_points(end);
    toc
end

kSpace_sb = apply_radial_ssg_traj_sp_SING(kSpace_mb, source_idx_all, kernel, nsms);
kSpace_sb = reshape(kSpace_sb, [sx, 210, 1, nc, nsms]);
N = NUFFT.init(kxt, kyt, 1, [6, 6], sx, sx);
test_im = NUFFT.NUFFT_adj(kSpace_sb, N);
test_im = squeeze(sos(test_im, 4));
figure
imagesc(test_im(:,:))
colormap gray
axis image
axis off
brighten(0.4)

[kSpace_cart(:, :, iframe, :, :), test] = apply_radial_ssg_traj_sp_SING(kSpace_mb, source_idx_all, kernel);

t1 = toc;
fprintf(repmat('\b',1,linelength));
fprintf(sprintf('calibration done, time = %.2f s\n', t1));
fprintf([repmat('-', [1, 75]), '\n'])

% display test image
if 0
    para.Recon.kSpace_size = [sx, sx];
    para.Recon.nor = nor_target;
    para.over_sampling = 1;
    para.Recon.sx = sx;
    
    w = ramp_filter_for_pre_interp(para);
    
    test_im = fftshift2(ifft2(fftshift2(temp_cart).*w));
    test_im = squeeze(sos(test_im, 4));
    
    figure
    imagesc(test_im(:,:))
    colormap gray
    axis image
    axis off
    brighten(0.4)
end