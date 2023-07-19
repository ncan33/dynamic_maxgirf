function img_all = section_kernel_apply_8_seg(img, kernel_im)

kSpace = fftshift2(fft2(fftshift2(img)));
img_all = zeros(size(img), 'like', img);

[sx, sy, ~] = size(img);
kernel_ex = 0;
%% 1
mask = ones(sx, sy);
mask(1:end,1:sy/2-kernel_ex,:,:,:) = 0;
mask(1:sx/2-kernel_ex,1:end,:,:,:) = 0;
[x, y] = meshgrid(1:sx*sqrt(2), 1:sy*sqrt(2));
x = imrotate(x - sx*sqrt(2)/2, 45);
x = crop_half_FOV(x,[sx,sy]);
x = x>0;
mask = mask.*x;
im_temp = fftshift2(ifft2(fftshift2(kSpace.*mask)));
im_temp = permute(im_temp,[1,2,3,4,6,5]) .* kernel_im{1};
im_temp = sum(im_temp, 4);
img_all = img_all + permute(im_temp, [1,2,3,5,6,4]);
%% 2
mask = ones(sx, sy);
mask(1:end,1:sy/2-kernel_ex,:,:,:) = 0;
mask(1:sx/2-kernel_ex,1:end,:,:,:) = 0;
[x, y] = meshgrid(1:sx*sqrt(2), 1:sy*sqrt(2));
x = imrotate(x - sx*sqrt(2)/2, 45);
x = crop_half_FOV(x,[sx,sy]);
x = x<0;
mask = mask.*x;
im_temp = fftshift2(ifft2(fftshift2(kSpace.*mask)));
im_temp = permute(im_temp,[1,2,3,4,6,5]) .* kernel_im{2};
im_temp = sum(im_temp, 4);
img_all = img_all + permute(im_temp, [1,2,3,5,6,4]);

%% 3
mask = ones(sx, sy);
mask(1:end,sy/2+1:end,:,:,:) = 0;
mask(1:sx/2,1:end,:,:,:) = 0;
[x, y] = meshgrid(1:sx*sqrt(2), 1:sy*sqrt(2));
x = imrotate(x - sx*sqrt(2)/2, -45);
x = crop_half_FOV(x,[sx,sy]);
x = x>0;
mask = mask.*x;
im_temp = fftshift2(ifft2(fftshift2(kSpace.*mask)));
im_temp = permute(im_temp,[1,2,3,4,6,5]) .* kernel_im{3};
im_temp = sum(im_temp, 4);
img_all = img_all + permute(im_temp, [1,2,3,5,6,4]);
%% 4
mask = ones(sx, sy);
mask(1:end,sy/2+1:end,:,:,:) = 0;
mask(1:sx/2,1:end,:,:,:) = 0;
[x, y] = meshgrid(1:sx*sqrt(2), 1:sy*sqrt(2));
x = imrotate(x - sx*sqrt(2)/2, -45);
x = crop_half_FOV(x,[sx,sy]);
x = x<0;
mask = mask.*x;
im_temp = fftshift2(ifft2(fftshift2(kSpace.*mask)));
im_temp = permute(im_temp,[1,2,3,4,6,5]) .* kernel_im{4};
im_temp = sum(im_temp, 4);
img_all = img_all + permute(im_temp, [1,2,3,5,6,4]);

%% 5
mask = ones(sx, sy);
mask(1:end,sy/2+1:end,:,:,:) = 0;
mask(sx/2+1:end,1:end,:,:,:) = 0;
[x, y] = meshgrid(1:sx*sqrt(2), 1:sy*sqrt(2));
x = imrotate(x - sx*sqrt(2)/2, 45);
x = crop_half_FOV(x,[sx,sy]);
x = x<0;
mask = mask.*x;
im_temp = fftshift2(ifft2(fftshift2(kSpace.*mask)));
im_temp = permute(im_temp,[1,2,3,4,6,5]) .* kernel_im{5};
im_temp = sum(im_temp, 4);
img_all = img_all + permute(im_temp, [1,2,3,5,6,4]);
%% 6
mask = ones(sx, sy);
mask(1:end,sy/2+1:end,:,:,:) = 0;
mask(sx/2+1:end,1:end,:,:,:) = 0;
[x, y] = meshgrid(1:sx*sqrt(2), 1:sy*sqrt(2));
x = imrotate(x - sx*sqrt(2)/2, 45);
x = crop_half_FOV(x,[sx,sy]);
x = x>0;
mask = mask.*x;
im_temp = fftshift2(ifft2(fftshift2(kSpace.*mask)));
im_temp = permute(im_temp,[1,2,3,4,6,5]) .* kernel_im{6};
im_temp = sum(im_temp, 4);
img_all = img_all + permute(im_temp, [1,2,3,5,6,4]);
%% 7
mask = ones(sx, sy);
mask(1:end,1:sy/2,:,:,:) = 0;
mask(sx/2+1:end,1:end,:,:,:) = 0;
[x, y] = meshgrid(1:sx*sqrt(2), 1:sy*sqrt(2));
x = imrotate(x - sx*sqrt(2)/2, -45);
x = crop_half_FOV(x,[sx,sy]);
x = x<0;
mask = mask.*x;
im_temp = fftshift2(ifft2(fftshift2(kSpace.*mask)));
im_temp = permute(im_temp,[1,2,3,4,6,5]) .* kernel_im{7};
im_temp = sum(im_temp, 4);
img_all = img_all + permute(im_temp, [1,2,3,5,6,4]);
%% 8
mask = ones(sx, sy);
mask(1:end,1:sy/2,:,:,:) = 0;
mask(sx/2+1:end,1:end,:,:,:) = 0;
[x, y] = meshgrid(1:sx*sqrt(2), 1:sy*sqrt(2));
x = imrotate(x - sx*sqrt(2)/2, -45);
x = crop_half_FOV(x,[sx,sy]);
x = x>0;
mask = mask.*x;
im_temp = fftshift2(ifft2(fftshift2(kSpace.*mask)));
im_temp = permute(im_temp,[1,2,3,4,6,5]) .* kernel_im{8};
im_temp = sum(im_temp, 4);
img_all = img_all + permute(im_temp, [1,2,3,5,6,4]);
