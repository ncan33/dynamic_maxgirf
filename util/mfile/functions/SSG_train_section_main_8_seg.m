function [kernel, kernel_im] = SSG_train_section_main_8_seg(ks, para)

[sx, sy, ~] = size(ks);
kernel_ex = 3;
%% 1
mask = ones(sx, sy);
mask(1:end,1:sy/2-kernel_ex,:,:,:) = 0;
mask(1:sx/2-kernel_ex,1:end,:,:,:) = 0;
[x, y] = meshgrid(1:sx*sqrt(2), 1:sy*sqrt(2));
x = imrotate(x - sx*sqrt(2)/2, 45);
x = crop_half_FOV(x,[sx,sy]);
x = x>0;
mask = mask.*x;
[kernel{1}, kernel_im{1}] = SSG_train_section(ks.*mask, para);
%% 2
mask = ones(sx, sy);
mask(1:end,1:sy/2-kernel_ex,:,:,:) = 0;
mask(1:sx/2-kernel_ex,1:end,:,:,:) = 0;
[x, y] = meshgrid(1:sx*sqrt(2), 1:sy*sqrt(2));
x = imrotate(x - sx*sqrt(2)/2, 45);
x = crop_half_FOV(x,[sx,sy]);
x = x<0;
mask = mask.*x;
[kernel{2}, kernel_im{2}] = SSG_train_section(ks.*mask, para);
%% 3
mask = ones(sx, sy);
mask(1:end,sy/2+1:end,:,:,:) = 0;
mask(1:sx/2,1:end,:,:,:) = 0;
[x, y] = meshgrid(1:sx*sqrt(2), 1:sy*sqrt(2));
x = imrotate(x - sx*sqrt(2)/2, -45);
x = crop_half_FOV(x,[sx,sy]);
x = x>0;
mask = mask.*x;
[kernel{3}, kernel_im{3}] = SSG_train_section(ks.*mask, para);
%% 4
mask = ones(sx, sy);
mask(1:end,sy/2+1:end,:,:,:) = 0;
mask(1:sx/2,1:end,:,:,:) = 0;
[x, y] = meshgrid(1:sx*sqrt(2), 1:sy*sqrt(2));
x = imrotate(x - sx*sqrt(2)/2, -45);
x = crop_half_FOV(x,[sx,sy]);
x = x<0;
mask = mask.*x;
[kernel{4}, kernel_im{4}] = SSG_train_section(ks.*mask, para);
%% 5
mask = ones(sx, sy);
mask(1:end,sy/2+1:end,:,:,:) = 0;
mask(sx/2+1:end,1:end,:,:,:) = 0;
[x, y] = meshgrid(1:sx*sqrt(2), 1:sy*sqrt(2));
x = imrotate(x - sx*sqrt(2)/2, 45);
x = crop_half_FOV(x,[sx,sy]);
x = x<0;
mask = mask.*x;
[kernel{5}, kernel_im{5}] = SSG_train_section(ks.*mask, para);
%% 6
mask = ones(sx, sy);
mask(1:end,sy/2+1:end,:,:,:) = 0;
mask(sx/2+1:end,1:end,:,:,:) = 0;
[x, y] = meshgrid(1:sx*sqrt(2), 1:sy*sqrt(2));
x = imrotate(x - sx*sqrt(2)/2, 45);
x = crop_half_FOV(x,[sx,sy]);
x = x>0;
mask = mask.*x;
[kernel{6}, kernel_im{6}] = SSG_train_section(ks.*mask, para);
%% 7
mask = ones(sx, sy);
mask(1:end,1:sy/2,:,:,:) = 0;
mask(sx/2+1:end,1:end,:,:,:) = 0;
[x, y] = meshgrid(1:sx*sqrt(2), 1:sy*sqrt(2));
x = imrotate(x - sx*sqrt(2)/2, -45);
x = crop_half_FOV(x,[sx,sy]);
x = x<0;
mask = mask.*x;
[kernel{7}, kernel_im{7}] = SSG_train_section(ks.*mask, para);
%% 8
mask = ones(sx, sy);
mask(1:end,1:sy/2,:,:,:) = 0;
mask(sx/2+1:end,1:end,:,:,:) = 0;
[x, y] = meshgrid(1:sx*sqrt(2), 1:sy*sqrt(2));
x = imrotate(x - sx*sqrt(2)/2, -45);
x = crop_half_FOV(x,[sx,sy]);
x = x>0;
mask = mask.*x;
[kernel{8}, kernel_im{8}] = SSG_train_section(ks.*mask, para);
