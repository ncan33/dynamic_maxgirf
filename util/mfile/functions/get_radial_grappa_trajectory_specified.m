function kernel = get_radial_grappa_trajectory_specified(Data, para)

[sx, nor_total, nc] = size(Data.kSpace);

% if isfield(para.Recon.ssg, 'calib_size')
%     calib_size = para.Recon.ssg.calib_size;
% else
%     calib_size = floor(sx/sqrt(2)/2)*2;
%     calib_size = [calib_size, calib_size];
%     para.Recon.ssg.calib_size = calib_size;
% end
% 
% if ~isfield(para.Recon.ssg, 'patch_size')
%     para.Recon.ssg.patch_size = [5, 5];
% end
% 
% if ~isfield(para.Recon.ssg, 'alpha')
%     para.Recon.ssg.patch_size = 4;
% end
% 
% if ~isfield(para.Recon.ssg, 'section')
%     para.Recon.ssg.section = 4;
% end

if ~isfield(para.Recon.rg, 'nor_calib')
    para.Recon.rg.nor_calib = 600;
end



theta = Data.theta;
[~, nor_target, nof_target] = size(Data.theta);
kx = reshape(Data.kx, [sx, nor_target, nof_target]);
ky = reshape(Data.ky, [sx, nor_target, nof_target]);

nor_calib = para.Recon.rg.nor_calib;

nof_total = floor(nor_total/nor_calib);
nor_total = nof_total * nor_calib;

Data.kSpace(:, nor_total+1:end, :, :, :) = [];
Data.kx(:, nor_total+1:end) = [];
Data.ky(:, nor_total+1:end) = [];

Data.kSpace = reshape(Data.kSpace, [sx, nor_calib, nof_total, nc]);
Data.kx = reshape(Data.kx, [sx, nor_calib, nof_total]);
Data.ky = reshape(Data.ky, [sx, nor_calib, nof_total]);

Data.N = NUFFT.init(Data.kx, Data.ky, 1, [6, 6], sx, sx);
Data.first_est = NUFFT.NUFFT_adj(Data.kSpace, Data.N);

%% set reconstruction parameters
% para.setting.ifplot = 1;
% para.setting.ifGPU = 0;
% para.Recon.epsilon = eps('single');
% para.Recon.step_size = 2;
% para.Recon.noi = 30;
% para.Recon.type = 'NUFFT coil';
% para.Recon.break = 1;
% para.Recon.weight_tTV = 0;
% para.Recon.weight_sTV = 0;

%% reconstruct aliasing-free coil images
% img = STCR_conjugate_gradient(Data, para);
% img = sliding_window_recon(Data);

kx_source = kx(:,:,1);
ky_source = ky(:,:,1);

[kx_target, ky_target] = meshgrid(-128:127, -128:127);

kSpace = reshape(Data.kSpace, [256, 72, 100, 9]);

kernel_size = [6, 6];
kSpace_out = train_radial_grappa_trajectory_specified(Data.first_est, kx_source, ky_source, kx_target, ky_target, kernel_size, Data.kSpace(:,1:72,1,:));

N = NUFFT.init(kx_target, ky_target, 1, [6, 6], sx, sx);
im = NUFFT.NUFFT_adj(kSpcae_out, N);


%% test
test_frame = 1:40;
kernel = train_radial_ssg_trajectory_specified(img, theta(:,:,test_frame), kx(:,:,test_frame), ky(:,:,test_frame), phase_mod_idx(:,:,test_frame));
kSpace = reshape(Data.kSpace, [sx, nor_target, nof_target, nc]);
kSpace_test = apply_radial_ssg_trajectroy_speficied(kSpace(:,:,test_frame,:), phase_mod_idx(:,:,test_frame), theta(:,:,test_frame), kernel);
N_test = NUFFT.init(kx(:,:,test_frame), ky(:,:,test_frame), 1, [6, 6], sx, sx);
image_test = NUFFT.NUFFT_adj(kSpace_test, N_test);
image_test = sum(image_test .* conj(Data.sens_map), 4);

phase_mod = reshape(Data.phase_mod, [1, nor_target, nof_target, 1, nsms]);

Data_test.kSpace = kSpace(:,:,test_frame,:);
Data_test.phase_mod_idx  = phase_mod_idx(:,:,test_frame);
Data_test.theta = theta(:,:,test_frame);
Data_test.first_est = image_test;
Data_test.sens_map = Data.sens_map;
Data_test.N = N_test;
Data_test.rssg = kernel;
Data_test.phase_mod = phase_mod(:,:,test_frame,:,:);

para.setting.ifplot = 1;
para.Recon.ssg.flag = 0;
para.Recon.rsg.flag = 1;
para.Recon.epsilon = eps('single');
para.Recon.step_size = 2;
para.Recon.noi = 150;
para.Recon.type = 'NUFFT RSG';
para.Recon.nSMS = 3;
para.Recon.break = 1;
scale = max(abs(Data.first_est(:)));
para.Recon.weight_tTV = scale * 0.02;
para.Recon.weight_sTV = scale * 0.001;

figure
imagesc([abs(image_test(:,:))*10; abs(recon(:,:))])
axis image
axis off
colormap gray
brighten(0.5)

