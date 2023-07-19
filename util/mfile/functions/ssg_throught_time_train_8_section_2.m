function [kernel, kernel_im] = ssg_throught_time_train_8_section_2(ks, para)

%% get mask
[sx, sy, ~] = size(ks);
kernel_ex = 3;

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
%% get kernel
for i=1:8
    [kernel{i}, kernel_im{i}] = SSG_tt_train_section(ks.*mask(:,:,i), para);
end

