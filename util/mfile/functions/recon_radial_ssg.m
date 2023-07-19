function recon_radial_ssg(para)

load(para.dir.file_name)

%% normalize k-space
kSpace = kSpace/max(abs(kSpace(:))) * 200^2*100;
nor_one = 24;
% para.setting.ifplot = 0;
% para.setting.ifGPU = 0;
para.Recon.nSMS = kSpace_info.usercv(14);
%%
[sx, nor_all, nc] = size(kSpace);

nof = floor(nor_all/nor_one);
nor_all = nof*nor_one;

kSpace(:,nor_all+1:end,:) = [];
kSpace = permute(kSpace,[1,2,4,3]);
kSpace_info.angle_mod(nor_all+1:end) = [];
kSpace_info.phase_mod(nor_all+1:end) = [];

theta = kSpace_info.angle_mod;
[kx ,ky] = get_k_coor(sx, theta, 0, sx/2+1);
%% compare RING 

% nor_temp = 300;
% N = NUFFT.init(kx(:,1:nor_temp), ky(:,1:nor_temp), 1, [6, 6], sx, sx);
% img_no_RING = NUFFT.NUFFT_adj(kSpace(:,1:nor_temp,:,:), N);
% 
correction = RING_SMS(kSpace, theta, para.Recon.nSMS);
para.Recon.RING = correction;
correction = correction.*permute([cos(theta);sin(theta)],[4,1,2,3]);
correction = squeeze(sum(correction,2));
kx = kx - correction(1,:,:);
ky = ky - correction(2,:,:);
% 
% N = NUFFT.init(kx(:,1:nor_temp), ky(:,1:nor_temp), 1, [6, 6], sx, sx);
% img_RING = NUFFT.NUFFT_adj(kSpace(:,1:nor_temp,:,:), N);
% 
% figure, imagesc([sos(img_no_RING), sos(img_RING), sos(img_no_RING) - sos(img_RING)]);
% colormap gray
% axis image
% brighten(0.6)
%% sens map
N = NUFFT.init(kx, ky, 1, [6,6], sx, sx);
phase_mod = kSpace_info.phase_mod;
phase_mod = exp(1i*2*pi*phase_mod/3);
phase_mod = cat(5, ones(size(phase_mod)), phase_mod, conj(phase_mod));
im = NUFFT.NUFFT_adj(kSpace .* conj(phase_mod), N);
Data.sens_map = get_sens_map(im, 'SMS');

%% train SSG kernel
Data.phase_mod = phase_mod;
Data.kSpace = kSpace;
Data.kx = kx;
Data.ky = ky;

para.Recon.ssg.nor_calib = 600;
para.Recon.ssg.calib_size = [140, 140];
para.Recon.ssg.patch_size = [5, 5];
para.Recon.ssg.alpha = 4;
para.Recon.ssg.section = 4;
% para.Recon.nSMS = 3;

[Data] = get_ssg_non_cart(Data, para);

%% recon multi frame
kx = reshape(kx, [sx, nor_one, nof]);
ky = reshape(ky, [sx, nor_one, nof]);
kSpace = reshape(kSpace, [sx, nor_one, nof, nc]);

%% get phase modulation
phase_mod = reshape(kSpace_info.phase_mod, [1, nor_one, nof]);
phase_mod = exp(1i*2*pi*phase_mod/3);
phase_mod = cat(5, ones(size(phase_mod)), phase_mod, conj(phase_mod));

%% get data structure 
nof_temp = min(300, nof);
nof_begin = 1;
nof_idx = nof_begin:nof_begin + nof_temp-1;

Data.kSpace = kSpace(:,:,nof_idx,:);
Data.phase_mod = phase_mod(:,:,nof_idx,:,:);
Data.N = NUFFT.init(kx(:,:,nof_idx,:), ky(:,:,nof_idx,:), 1, [6, 6], sx, sx);
Data.first_est = NUFFT.NUFFT_adj(Data.kSpace .* conj(Data.phase_mod), Data.N)/para.Recon.nSMS;
Data.first_est = section_kernel_apply(Data.first_est, Data.ssg);
% Data.sens_map = get_sens_map(Data.first_est, 'SMS');
Data.first_est = sum(Data.first_est .* conj(Data.sens_map), 4);
% Data.sens_map_conj = sens_conj_ssg;
% Data.ssg = kernel_im;
% Data.kernel = sens_conj_ssg;

%% get para
para.Recon.ssg = 1;

% para.setting.ifplot = 0;
% para.setting.ifGPU = 1;
para.Recon.epsilon = eps('single');
para.Recon.step_size = 2;
para.Recon.noi = 150;
para.Recon.type = 'NUFFT';
% para.Recon.nSMS = 3;
para.Recon.break = 1;
scale = max(abs(Data.first_est(:)));
para.Recon.weight_tTV = scale * 0.3;
para.Recon.weight_sTV = scale * 0.002;
%% recon with ssg
clearvars -except Data para file_name
% img = STCR_gradient_descent(Data, para);
[Image, para] = STCR_conjugate_gradient(Data, para);
Image = abs(Image);

calib_size = para.Recon.ssg.calib_size;
patch_size = para.Recon.ssg.patch_size;
folder = strfind(file_name, '/RawData/');
folder = file_name(1:folder);
P = file_name(folder+9:end-4);
if isempty(dir(sprintf('%sReconData', folder)))
    mkdir(sprintf('%sReconData', folder))
end
file_name = sprintf('%sReconData/%s_%g_ray_calibration_size_%g_%g_segments_%g_a_%g_kernel_size_%g_%g.mat', folder, P, para.Recon.ssg.nor_calib, calib_size(1), calib_size(2), para.Recon.ssg.section, para.Recon.ssg.alpha, patch_size(1), patch_size(2));
save(file_name,'Image', 'para');

%% recon without ssg
para.Recon.ssg = 0;
[Image, para] = STCR_conjugate_gradient(Data, para);
Image = abs(Image);
file_name = sprintf('%sReconData/%s.mat', folder, P);
save(file_name,'Image', 'para');
