function [kernel, kernel_im_space] = SSG_train(ks, para)

% ks is single-band calibration data
[sx, sy, nc, nsms, nslice] = size(ks);
ps = para.patch_size;

kt = zeros([sx, sy, nc, nslice], 'like', ks);
nx = sx - ps(1) + 1;
ny = sy - ps(2) + 1;
npatch = nx*ny;

% sum single slice with aliased SMS slices data
% target is single slice data
ks_sum_sms = sum(ks, 4);
for i=1:nsms
    kt(:,:,:,i) = ks(:,:,:,i,i);
    ks(:,:,:,i,i) = ks_sum_sms(:,:,:,i);
end

kernel = zeros([nc, prod(ps)*nc, nsms], 'like', ks);
for i=1:nslice
    source = zeros([prod(ps)*nc, npatch*nsms], 'like', ks);
    target = zeros([nc, npatch*nsms], 'like', ks);
    for j=1:nsms
        idx = (1:npatch) + (j-1)*npatch;
        source_temp = get_patch_sl1(ks(:,:,:,j,i), ps);
        % the training aims to noun the aliaed k-space data, so the source
        % is zero 
        if i==j
            target_temp = get_patch_center_sl1(kt(:,:,:,i), ps);
        else
            target_temp = zeros([nc, npatch], 'like', ks);
        end
        source(:,idx) = source_temp;
        target(:,idx) = target_temp;
    end
    kernel(:,:,i) = target * pinv(source);
end

if nargout > 1
    kernel = reshape(kernel, [nc, ps, nc, nsms]);
    kernel = permute(kernel, [2, 3, 6, 4, 1, 5]);
    kernel = rot90(kernel, 2);
    
    kernel_im_space = zeros([para.imsize, 1, nc, nc, nsms]);
    idx_x = round((para.imsize(1) - ps(1))/2);
    idx_x = idx_x + 1 : idx_x + ps(1);
    idx_y = round((para.imsize(2) - ps(2))/2);
    idx_y = idx_y + 1 : idx_y + ps(2);
    kernel_im_space(idx_x, idx_y, :,:,:,:) = kernel;
    
    kernel_im_space = ifftshift2(kernel_im_space);
    kernel_im_space = ifft2(kernel_im_space);
    kernel_im_space = fftshift2(kernel_im_space);
%     kernel_im_space = rot90(kernel_im_space, 2);
    kernel_im_space = kernel_im_space * prod(para.imsize);
end