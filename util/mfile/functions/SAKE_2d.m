function im = SAKE_2d(kspace, N, kernel_size, wnRank, nIter, show)
%
% res = SAKE(DATA, kSize, wnRank,nIter)
%
% function performs an autocalibration reconstrution without calibration
% lines based on the SAKE method (P. Shin et. al, "Calibrationless Parallel 
% Imaging Reconstruction Based on Structured Low-Rank Matrix Completion"
% 2013, submitted to MRM. 
%
% It is recommended that this method be used to complete a calibration
% region and then use ESPIRiT to generate ESPIRiT maps. 
%
% INPUTS:
%       DATA [nx X ny X nc]  - 2D multi-coil k-space data with missing
%                               entries, preferably random! missing entries
%                               should be EXACTLY = 0!
%       kSize [ 1 x 2]       - sliding window size
%       wnRank [ 1 x 1]       - the window-normalized rank to enforce 
%                               (# of coils = full rank, 2.0 would be typical 
%                               for 8 coils with window size [5 5]
%       nIter                - number of iTerations.
%       show                 - show intermediate reconstruction show=0 will
%                               skip. Show = 100 will plot in figure 100
%
%
% Outputs:
%       res -                - 2D multi coil k-space data where the missing
%                               data was filled. 
%
% (c) Michael Lustig 2012
%



% mask = abs(kspace)>0;

% res = kspace;
kspace          = single(kspace);
im              = NUFFT.NUFFT_adj(kspace, N);
kspace_cart     = fftshift2(fft2(fftshift2(im)));
[sx, sy, nc]    = size(kspace_cart);
keep            = 1:floor(wnRank * prod(kernel_size));
mask            = circular_mask(sx, sx/2);
us_mask         = logical(abs(kspace(1, :, :, 1)));

kspace_cart = kspace_cart .* mask;

kspace_cart = gpuArray(kspace_cart);
N.S         = gpuArray(N.S);
kspace      = gpuArray(kspace);
 tic
 figure(1)
    clf
    imagesc(abs(sos(im)))
    axis image
    drawnow
    
for n = 1:nIter
   
    % reorder data to get Hankel structure. 
    tmp = single(im2row(squeeze(kspace_cart), kernel_size)); 
    [tsx, tsy, tsz] = size(tmp);
    A = reshape(tmp, tsx, tsy*tsz);
    
    % SVD thresholding
    [U, S, V] = svdecon(A);
    A = U(:, keep) * S(keep,keep) * V(:,keep)';
    
    % Enforce Hankel structure
    A = reshape(A, tsx, tsy, tsz);
    tmp = single(row2im(A, [sx, sy, nc], kernel_size));
    
    % enforce data consistency
    im = permute(fftshift2(ifft2(fftshift2(tmp))), [1, 2, 4, 3]);
    
    f_update = NUFFT.NUFFT_adj(kspace - NUFFT.NUFFT(im, N) .* us_mask, N);
    
    im = im + f_update * 0.5;
    kspace_cart = fftshift2(fft2(fftshift2(im))) .* mask;
    
%     res = tmp.*(1-mask) + kspace;

if mod(n, 30) == 0
    toc
    tic
    figure(1)
    clf
    imagesc(abs(sos(im)))
    axis image
    drawnow
end
end

%{
for n = 1:nIter
    tic
    % reorder data to get Hankel structure. 
    tmp = im2row(squeeze(kspace_cart), kernel_size); 
    [tsx, tsy, tsz] = size(tmp);
    A = reshape(tmp, tsx, tsy*tsz);
    
    % SVD thresholding
    [U, S, V] = svd(A, 'econ');
    
    A = U(:, keep) * S(keep,keep) * V(:,keep)';
    
    % Enforce Hankel structure
    A = reshape(A, tsx, tsy, tsz);
    tmp = row2im(A, [sx, sy, nc], kernel_size);
    
    % enforce data consistency
    s_update = permute(fftshift2(ifft2(fftshift2(tmp))), [1, 2, 4, 3]) - im;
    
    f_update = NUFFT.NUFFT_adj(kspace - NUFFT.NUFFT(im, N), N);
    
    im = im + f_update * 0.3 + s_update;
    kspace_cart = fftshift2(fft2(fftshift2(im))) .* mask;
    
%     res = tmp.*(1-mask) + kspace;
    figure(1)
    imagesc(sos(im))
    axis image
    drawnow
    toc
end

%}

