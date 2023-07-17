%% add ismrmd and mfile path
addpath /server/home/ncan/ismrmrd
addpath ./matlab/

ccc
if ~isfolder('./recon_data')
    mkdir ./recon_data
end

all_dat     = dir('/server/sdata/ncan/mri_data/disc/lung/vol0457_20221021/raw_hawk/usc_disc_yt_2022_10_21_133643_dual-te_dynamic.mat');
nfile       = length(all_dat);
% naverage    = 1;
narm_frame  = 2;
para.Recon.narm = narm_frame;   % number of arms per frame
para.Recon.time_frames = 'all'; % set to 'all' for reconstructructing all frames. specify certain range if you wish, i.e., 1:100
% nprep       = 5;
ifsave      = 1;
file_index = 1;

%% read data
file_name = fullfile(all_dat(file_index).folder, all_dat(file_index).name);
fprintf('--------------------------------- \n')
fprintf('reading data:  %d/%d \n', file_index, nfile)
fprintf('file name:     %s \n', file_name)
fprintf('\n')

load(file_name)

narm    = size(kspace, 2);
nview   = kspace_info.user_interleaves;

%% reconstruction
tic

view_order = kspace_info.viewOrder;

echo_idx = view_order > nview;
echo_idx = echo_idx + 1;

%% size
nsample         = size(kspace, 1);
ncoil           = size(kspace, 3);

%% split kspace
kspace_echo_1 = kspace(:, echo_idx == 1, :);
kspace_echo_2 = kspace(:, echo_idx == 2, :);
%clear kspace

narm_total      = min(size(kspace_echo_1, 2), size(kspace_echo_2, 2));

%% orgnize the data to frames
nframes = floor(narm_total / narm_frame);

narm_total = nframes * narm_frame;

kspace_echo_1(:, narm_total + 1 : end, :) = [];
kspace_echo_2(:, narm_total + 1 : end, :) = [];

kspace_echo_1  = reshape(kspace_echo_1, [nsample, narm_frame, nframes, ncoil]);
kspace_echo_2  = reshape(kspace_echo_2, [nsample, narm_frame, nframes, ncoil]);

%% kspace trajectory and view order
view_order_echo_1 = view_order(echo_idx == 1);
view_order_echo_2 = view_order(echo_idx == 2);
view_order_echo_1(narm_total + 1 : end) = [];
view_order_echo_2(narm_total + 1 : end) = [];
view_order_echo_1  = reshape(view_order_echo_1, [narm_frame, nframes]);
view_order_echo_2  = reshape(view_order_echo_2, [narm_frame, nframes]);

kx = kspace_info.kx_GIRF;
ky = kspace_info.ky_GIRF;

kx_echo_1 = zeros(size(kspace_echo_1,1), size(kspace_echo_1,2), size(kspace_echo_1,3));
ky_echo_1 = kx_echo_1;
%kx_echo_1 = zeros(nsample, narm_frame, nframes);
kx_echo_2 = zeros(size(kspace_echo_2,1), size(kspace_echo_2,2), size(kspace_echo_2,3));
ky_echo_2 = kx_echo_2;
%kx_echo_1 = zeros(nsample, narm_frame, nframes);

for i = 1:narm_frame
    for j = 1:nframes
        kx_echo_1(:,i,j) = kx(:, view_order_echo_1(i,j));
        kx_echo_2(:,i,j) = kx(:, view_order_echo_2(i,j) - 10);
        
        ky_echo_1(:,i,j) = ky(:, view_order_echo_1(i,j));
        ky_echo_2(:,i,j) = ky(:, view_order_echo_2(i,j) - 10);
    end
end

%% matrix size
res = [kspace_info.user_ResolutionX, kspace_info.user_ResolutionY];
fov = [kspace_info.user_FieldOfViewX, kspace_info.user_FieldOfViewY];
matrix_size_keep    = round(fov ./ res / 2) * 2;
matrix_size         = round(matrix_size_keep * 2 / 2) * 1.5;
para.Recon.matrix_size = matrix_size;

%% select time frames - echo 1
if isnumeric(para.Recon.time_frames)
    time_frames = para.Recon.time_frames;
    
    kx_echo_1 = kx_echo_1(:, :, time_frames);
    ky_echo_1 = ky_echo_1(:, :, time_frames);
    kspace_echo_1 = kspace_echo_1(:, :, time_frames, :);
    
    kx_echo_2 = kx_echo_2(:, :, time_frames);
    ky_echo_2 = ky_echo_2(:, :, time_frames);
    kspace_echo_2 = kspace_echo_2(:, :, time_frames, :);
end

para.Recon.FOV = 4.8; %units of decimeter for some reason

Data.N = NUFFT.init(kx_echo_1*para.Recon.FOV, ky_echo_1*para.Recon.FOV, 1, [4, 4], para.Recon.matrix_size(1)*para.Recon.FOV, para.Recon.matrix_size(1)*para.Recon.FOV);
Data.N.W = kspace_info.DCF(:, 1);

Data.kSpace = kspace_echo_1;
Data.first_est = NUFFT.NUFFT_adj(Data.kSpace, Data.N);

scale = max(abs(Data.first_est(:)));

Data.sens_map = get_sens_map(Data.first_est, '2D');
Data.first_est = sum(Data.first_est .* conj(Data.sens_map), 4);

%% set parameters
para.weight_tTV = 0; % CHANGE THIS CHANGE THIS CHANGE THIS CHANGE THIS
para.weight_sTV = 0; % CHANGE THIS CHANGE THIS CHANGE THIS CHANGE THIS
para.Recon.weight_tTV = scale * para.weight_tTV; % temporal regularization weight
para.Recon.weight_sTV = scale * para.weight_sTV; % spatial regularization weight

para.setting.ifplot = 1;        % display image and cost during reconstruction
para.setting.ifGPU = 1;         % set to 1 when you want to use GPU

para.Recon.no_comp = ncoil;     % number of coils
%para.Recon.FOV = 1.25;         % reconstruction FOV

para.Recon.epsilon = eps('single');
para.Recon.step_size = 2;
para.Recon.ifContinue = 0;
para.Recon.noi = 150; % number of iterations
para.Recon.type = '2D Spiral server';
para.Recon.break = 1;

%clearvars -except Data para

%% conjugate gradient reconstruction - echo 1
[Image_recon, para] = STCR_conjugate_gradient(Data, para);
Image_recon = rot90(Image_recon);
im_echo_1 = crop_half_FOV(Image_recon, para.Recon.matrix_size);

%% select time frames - echo 2
Data.kSpace = kspace_echo_2;
Data.N = NUFFT.init(kx_echo_1*para.Recon.FOV, ky_echo_1*para.Recon.FOV, 1, [4, 4], para.Recon.matrix_size(1)*para.Recon.FOV, para.Recon.matrix_size(1)*para.Recon.FOV);

Data.first_est = NUFFT.NUFFT_adj(Data.kSpace, Data.N);

scale = max(abs(Data.first_est(:)));

Data.sens_map = get_sens_map(Data.first_est, '2D');
Data.first_est = sum(Data.first_est .* conj(Data.sens_map), 4);

%% conjugate gradient reconstruction - echo 2
[Image_recon, para] = STCR_conjugate_gradient(Data, para);
Image_recon = rot90(Image_recon);
im_echo_2 = crop_half_FOV(Image_recon, para.Recon.matrix_size);

%% end excerpt

%clearvars -except im_echo_1 im_echo_2 kspace_info

%{
%%save 
f_mag = imageMRI(sos([im_echo_1(:, :, end, :), im_echo_2(:, :, end, :)]));
f_phase = imageMRI(angle([im_echo_1(:, :, end, 1), im_echo_2(:, :, end, 1)]));
if ifsave
    save_name   = sprintf('./recon_data/%s_recon.mat', all_dat(file_index).name(1:end-8));
%         figure_name = sprintf('./figure/%s_recon_naverage_%g.png', all_dat(i).name(1:end-8), naverage);
%         saveas(f, figure_name);
    save(save_name, 'im_echo_1', 'im_echo_2', 'kspace_info', '-v7.3');
end

toc

%}