function kSpace_out = train_radial_grappa_trajectory_specified(im_sb, kx_source, ky_source, kx_target, ky_target, kernel_size, kSpace_in)

[sx, sy, nof, nc] = size(im_sb);
kSpace = fft2(im_sb);
target_size = size(kx_target);
source_size = size(kx_source);

kx_target = kx_target(:);
ky_target = ky_target(:);

kx_source = kx_source(:);
ky_source = ky_source(:);

n_target = length(kx_target);

kSpace_in = reshape(kSpace_in, [prod(source_size), nc]);

kSpace_out = zeros([prod(target_size), nc], 'single');

parfor ii = 1:n_target
    fprintf(sprintf('%g/%g\n', ii, n_target))
    kx_t = kx_target(ii);
    ky_t = ky_target(ii);
    
    % find all source points within kernel size
    idx = sos([kx_t - kx_source, ky_t - ky_source]) < kernel_size(1);
    np = sum(idx);
    
    if np >= 1
        kx_s = kx_source(idx);
        ky_s = ky_source(idx);
        
        kx_ = [kx_t; kx_s];
        ky_ = [ky_t; ky_s];
        
        N = NUFFT.init(kx_, ky_, 1, [6, 6], sx, sx);
        
        kSpace_calib = NUFFT.cart2rad(permute(kSpace, [1, 2, 5, 3, 4]), N);
        
        target = kSpace_calib(1, :, :);
        target = reshape(target, [nc, nof]);
        
        source = kSpace_calib(2:end, :, :);
        source = reshape(source, [np*nc, nof]);
        
        kernel = target * pinv(source);
        kSpace_source = kSpace_in(idx, :);
        kSpace_source = kSpace_source(:);
        kSpace_target = kernel * kSpace_source;
        kSpace_out(ii, :) = kSpace_target;
    end
end
kSpace_out = reshape(kSpace_out, [target_size, 1, nc]);

keyboard