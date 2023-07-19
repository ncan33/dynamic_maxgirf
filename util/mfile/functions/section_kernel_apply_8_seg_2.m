function img_all = section_kernel_apply_8_seg_2(img, kernel_im)


%% get mask
[sx, sy, ~] = size(img);
kernel_ex = 0;

mask_0 = zeros(sx, sy, 'single');
mask_0(1:sx/2, 1:sy/2) = 1;
mask_0 = circshift(mask_0, [sx/4, sy/4]);

mask_1 = ones(sx, sy, 'single');
mask_1(1:end,1:sy/2-kernel_ex) = 0;
mask_1(1:sx/2-kernel_ex,1:end) = 0;
mask(:,:,1) = mask_1.*mask_0;

mask_1 = ones(sx, sy, 'single');
mask_1(1:end,1:sy/2-kernel_ex) = 0;
mask_1(1:sx/2-kernel_ex,1:end) = 0;
mask(:,:,2) = mask_1.*~mask_0;

mask_1 = ones(sx, sy, 'single');
mask_1(1:end,sy/2+1:end) = 0;
mask_1(1:sx/2,1:end) = 0;
mask(:,:,3) = mask_1.*mask_0;

mask_1 = ones(sx, sy, 'single');
mask_1(1:end,sy/2+1:end) = 0;
mask_1(1:sx/2,1:end) = 0;
mask(:,:,4) = mask_1.*~mask_0;

mask_1 = ones(sx, sy, 'single');
mask_1(1:end,sy/2+1:end) = 0;
mask_1(sx/2+1:end,1:end) = 0;
mask(:,:,5) = mask_1.*mask_0;

mask_1 = ones(sx, sy, 'single');
mask_1(1:end,sy/2+1:end) = 0;
mask_1(sx/2+1:end,1:end) = 0;
mask(:,:,6) = mask_1.*~mask_0;

mask_1 = ones(sx, sy, 'single');
mask_1(1:end,1:sy/2) = 0;
mask_1(sx/2+1:end,1:end) = 0;
mask(:,:,7) = mask_1.*mask_0;

mask_1 = ones(sx, sy, 'single');
mask_1(1:end,1:sy/2) = 0;
mask_1(sx/2+1:end,1:end) = 0;
mask(:,:,8) = mask_1.*~mask_0;

%% apply kernel
kSpace = fftshift2(fft2(fftshift2(img)));
img_all = zeros(size(img), 'like', img);

for i=1:8
    im_temp = fftshift2(ifft2(fftshift2(kSpace.*mask(:,:,i))));
    im_temp = permute(im_temp,[1,2,3,4,6,5]) .* kernel_im{i};
    im_temp = sum(im_temp, 4);
    
%     im_temp = fftshift2(fft2(fftshift2(im_temp)));
%     im_temp = im_temp.*mask(:,:,i);
%     im_temp = fftshift2(ifft2(fftshift2(im_temp)));
    
    img_all = img_all + permute(im_temp, [1,2,3,5,6,4]);
end

