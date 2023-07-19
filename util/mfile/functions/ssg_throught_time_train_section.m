function [kernel, kernel_im] = ssg_throught_time_train_section(ks, para)

[sx, sy, ~] = size(ks);
kernel_ex = 0;

ks_temp = ks;
ks_temp(1:end,1:sy/2-kernel_ex,:,:,:,:) = 0;
ks_temp(1:sx/2-kernel_ex,1:end,:,:,:,:) = 0;
[kernel{1}, kernel_im{1}] = SSG_tt_train_section(ks_temp, para);

ks_temp = ks;
ks_temp(1:end,sy/2+1+kernel_ex:end,:,:,:,:) = 0;
ks_temp(1:sx/2-kernel_ex,1:end,:,:,:,:) = 0;
[kernel{2}, kernel_im{2}] = SSG_tt_train_section(ks_temp, para);

ks_temp = ks;
ks_temp(1:end,sy/2+1+kernel_ex:end,:,:,:,:) = 0;
ks_temp(sx/2+1+kernel_ex:end,1:end,:,:,:,:) = 0;
[kernel{3}, kernel_im{3}] = SSG_tt_train_section(ks_temp, para);

ks_temp = ks;
ks_temp(1:end,1:sy/2-kernel_ex,:,:,:,:) = 0;
ks_temp(sx/2+1+kernel_ex:end,1:end,:,:,:,:) = 0;
[kernel{4}, kernel_im{4}] = SSG_tt_train_section(ks_temp, para);