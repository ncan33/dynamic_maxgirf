function [kspace_out, kernel] = GRAPPA_train(kspace, para)

% kspace is all k-space data
[sx, sy, nc] = size(kspace);
ps = para.patch_size;

% kt is kspace target
% kt = zeros([sx, sy, nslice], 'like', kspace);

sy_calib = sum(para.mask_calib(1, :));
kspace_calib = kspace(repmat(para.mask_calib, [1, 1, nc]));
kspace_calib = reshape(kspace_calib, [sx, sy_calib, nc]);

nx = sx - ps(1) + 1;
ny = sy_calib - ps(2) + 1;
npatch = nx * ny;


kernel = zeros([nc, prod(ps)*nc], 'like', kspace);

mask_undersample = para.mask(1:ps(1), 1:ps(2));
source = get_patch_with_mask(kspace_calib, ps, mask_undersample);
target = get_patch_center_sl1(kspace_calib, ps);

kernel = target * pinv(source);


source = get_patch_with_mask(kspace, ps, true(3));
target = kernel * source;
kspace_temp = reshape(target, [nc, sx - (ps(1) - 1), sy - (ps(2) - 1)]);
kspace_temp = permute(kspace_temp, [3, 2, 1]);
kspace_temp(end + ps(1) - 1, end + ps(2) - 1, :) = 0;
kspace_temp = circshift(kspace_temp, [(ps(1) - 1)/2, (ps(2) - 1)/2, 0]);

kspace_out = kspace + kspace_temp .* ~para.mask;

%% kernel in image space

if nargout > 1
    kernel = reshape(kernel, [nc, ps, nc]);
    kernel = permute(kernel, [2, 3, 4, 1]);
    kernel = rot90(kernel, 2);
    
    kernel_im_space = zeros([sx, sx, nc, nc]);
    idx_x = round((sx - ps(1))/2);
    idx_x = idx_x + 1 : idx_x + ps(1);
    idx_y = round((sy - ps(2))/2);
    idx_y = idx_y + 1 : idx_y + ps(2);
    kernel_im_space(idx_x, idx_y, :,:,:,:) = kernel;
    
    kernel_im_space = ifftshift2(kernel_im_space);
    kernel_im_space = ifft2(kernel_im_space);
    kernel_im_space = fftshift2(kernel_im_space);
%     kernel_im_space = rot90(kernel_im_space, 2);
    kernel_im_space = kernel_im_space * sx * sy;
end


%% 