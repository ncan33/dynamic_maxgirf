function [img_all, img_only_apply_k_in_center] = section_kernel_apply(img, kernel_im)

%% single section
if size(kernel_im, 1) > 1
    img_all = permute(img, [1,2,3,4,6,5]) .* kernel_im;
    img_all = sum(img_all, 4);
    img_all = permute(img_all, [1,2,3,5,6,4]);
    return
end

%% get mask
[sx, sy, ~] = size(img);

mask = ones(sx, sy, 4, 'single');
mask(1:end,1:sy/2,1) = 0;
mask(1:sx/2,1:end,1) = 0;

mask(1:end,sy/2+1:end,2) = 0;
mask(1:sx/2,1:end,2) = 0;

mask(1:end,sy/2+1:end,3) = 0;
mask(sx/2+1:end,1:end,3) = 0;

mask(1:end,1:sy/2,4) = 0;
mask(sx/2+1:end,1:end,4) = 0;

%% apply kernel
kSpace = fftshift2(fft2(fftshift2(img)));
img_all = zeros(size(img), 'like', img);
for i=1:4
    im_temp = fftshift2(ifft2(fftshift2(kSpace.*mask(:,:,i)))); clear k_temp
    im_temp = permute(im_temp,[1,2,3,4,6,5]) .* kernel_im{i};
    im_temp = sum(im_temp, 4);
    img_all = img_all + permute(im_temp, [1,2,3,5,6,4]);
end

%% 
if 0
    mask_center = zeros([sx, sy], 'single');
    mask_center(sx/2, sy/2) = 1;
    mask_center = bwdist(mask_center);
    mask_center(mask_center>sx/8) = 0;
    mask_center(sx/2, sy/2) = 1;
    mask_center = logical(mask_center);
    
    k_after_k = fftshift2(fft2(fftshift2(img_all)));
    k_only_apply_k_in_center = k_after_k.*mask_center + kSpace.*~mask_center;
    img_only_apply_k_in_center = fftshift2(ifft2(fftshift2(k_only_apply_k_in_center)));
    
end
