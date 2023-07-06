clear all; close all; clc;

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
narm_frame  = 10;
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
clear kspace

narm_total      = min(size(kspace_echo_1, 2), size(kspace_echo_2, 2));

%% orgnize the data to frames
nframes = floor(narm_total / narm_frame);

narm_total = nframes * narm_frame;

kspace_echo_1(:, narm_total + 1 : end, :) = [];
kspace_echo_2(:, narm_total + 1 : end, :) = [];

kspace_echo_1  = reshape(kspace_echo_1, [nsample, narm_frame, nframes, ncoil]);
kspace_echo_2  = reshape(kspace_echo_2, [nsample, narm_frame, nframes, ncoil]);

%% begin recon_stcr_VD_spiral
addpath /server/home/ytian/mfile/
addpath /server/home/ytian/mfile/functions/
addpath /server/home/ytian/mfile/registrtation/
addpath /server/home/ytian/mfile/quantification/
addpath /server/home/ytian/mfile/ModelBasedMy/
addpath /server/home/ytian/mfile/mapVBVD/
addpath /server/home/ytian/mfile/vdspiral/
addpath ./

FOV_recon = [480, 480]; % mm
% Narms_per_frame = 3;
% recon_arms = (1:89*2) + 23; %'all'; %100:1647; %100:279;
% temporal_resolution = 500; % ms

%% set settings
% para.setting.off_res_corr   = 1;
para.setting.ifplot         = 1;        % plot convergence during reconstruction 
para.setting.ifGPU          = 1;        % set to 1 when you want to use GPU

%% set recon parameters
para.weight_tTV = 0; % JT's is 0.04;      % temporal TV regularizaiton parameter (normalized by F^T d)
para.weight_sTV = 0;                        % spatial TV regularizaiton parameter (normalized by F^T d)

para.Recon.epsilon      = eps('single');        % small vale to avoid singularity in TV constraint
para.Recon.step_size    = 2;                    % initial step size
para.Recon.noi          = 20;                  % number of CG iterations
para.Recon.type         = '2D Spiral server';   % 2D spiral
para.Recon.break        = 1;                    % stop iteration if creteria met. Otherwise will run to noi

% parameter for low-rank applying off-res map
% para.L              = 6;
% para.n_hist_bins    = 40;

%% echo 1

% para.te                 = 0.592 / 1000; % echo time
% para.crop_thresh        = 1 / para.te / 2;
% para.NarmsFull          = kspace_info.user_interleaves;
% para.sizehanningfilter  = 21;
% para.tdwell             = 1 ./ kspace_info.user_samplingRate / 1000; % dwell time
% para.nsample            = kspace_info.user_samples;
% para.tread              = para.nsample * para.tdwell;
% para.t                  = para.te + [0:para.nsample-1] * para.tread / para.nsample;

TR = kspace_info.user_TR / 1000;
temporal_resolution = TR*20; % units of ms. For the vol0457 data, TR*4 = 96.24 ms
para.Recon.narm = floor(temporal_resolution / TR);
Narms_per_frame = para.Recon.narm;

res = [kspace_info.user_ResolutionX, kspace_info.user_ResolutionY];

para.Recon.TR  = TR;
para.Recon.res = res;
para.kspace_info = kspace_info;

kspace = permute(kspace_echo_1, [1, 2, 4, 3]);

matrix_size = round(FOV_recon ./ res / 2) * 2;

para.Recon.image_size = matrix_size;

matrix_size_keep = [kspace_info.user_FieldOfViewX, kspace_info.user_FieldOfViewX] ./ res;
para.Recon.matrix_size_keep = round(matrix_size_keep);


kx = kspace_info.kx_GIRF * matrix_size(1);
ky = kspace_info.ky_GIRF * matrix_size(2);
    
% kx = kspace_info.kx * matrix_size(1);
% ky = kspace_info.ky * matrix_size(2);
    
% Nsamples = kspace_info.user_samples * 2;

%% gradient delay correction
% kx = kx(2:end, :);
% ky = ky(2:end, :);
% kspace = kspace(1:end-1, :, :, :);
% kspace_info.DCF = kspace_info.DCF(2:end, :);
    
viewOrder = kspace_info.viewOrder;

if exist('recon_arms', 'var')
    if isnumeric(recon_arms)
        kspace = kspace(:, recon_arms, :, :);
        viewOrder = viewOrder(recon_arms);
    else
        kspace = kspace(:, 151:end, :, :);
        viewOrder = viewOrder(151:end);
    end
end

GA_steps = size(kx, 2);
Narms_total = size(kspace, 2);
Nframes = floor(Narms_total / Narms_per_frame);
Narms_total = Nframes * Narms_per_frame;
Ncoil = size(kspace, 4);
Nsample = size(kspace, 1);

kx = repmat(kx, [1, ceil(Narms_total / GA_steps)]);
ky = repmat(ky, [1, ceil(Narms_total / GA_steps)]);

kspace(:, Narms_total + 1 : end, :, :) = [];
viewOrder(Narms_total + 1 : end) = [];

kx = kx(:, viewOrder);
ky = ky(:, viewOrder);

kspace = reshape(kspace, [Nsample, Narms_per_frame, Nframes, Ncoil]);

Nsample_k = size(kx, 1);
kx = reshape(kx, [Nsample_k, Narms_per_frame, Nframes]);
ky = reshape(ky, [Nsample_k, Narms_per_frame, Nframes]);
    
%% DCF
% kspace_info.DCF(kspace_info.DCF > 0.2) = 0.2;
% kspace_info.DCF = kspace_info.DCF * 5;
    
%% spiral out
delay = 0;

kx_in = kx((1:floor(Nsample/2)) + delay, :, :);
ky_in = ky((1:floor(Nsample/2)) + delay, :, :);

DCF = kspace_info.DCF((1:floor(Nsample/2)) + delay, :);

Data.N = NUFFT.init(kx_in, ky_in, 1, [6, 6], matrix_size(1), matrix_size(1));
Data.N.W = DCF(:, 1);

Data.kSpace = kspace(1:floor(Nsample/2), :, :, :);
Data.first_est = NUFFT.NUFFT_adj(Data.kSpace, Data.N);
Data.sens_map = get_sens_map(Data.first_est, '2D');
Data.first_est = sum(Data.first_est .* conj(Data.sens_map), 4);

scale = max(abs(Data.first_est(:)));

para.Recon.no_comp = Ncoil;
para.Recon.weight_tTV = scale * para.weight_tTV; % temporal regularization weight
para.Recon.weight_sTV = scale * para.weight_sTV; % spatial regularization weight

[Image_1, para] = STCR_conjugate_gradient(Data, para);

%% echo 2
kspace = permute(kspace_echo_2, [1, 2, 4, 3]);

matrix_size = round(FOV_recon ./ res / 2) * 2;

para.Recon.image_size = matrix_size;

matrix_size_keep = [kspace_info.user_FieldOfViewX, kspace_info.user_FieldOfViewX] ./ res;
para.Recon.matrix_size_keep = round(matrix_size_keep);


kx = kspace_info.kx_GIRF * matrix_size(1);
ky = kspace_info.ky_GIRF * matrix_size(2);
    
% kx = kspace_info.kx * matrix_size(1);
% ky = kspace_info.ky * matrix_size(2);
    
% Nsamples = kspace_info.user_samples * 2;

%% gradient delay correction
% kx = kx(2:end, :);
% ky = ky(2:end, :);
% kspace = kspace(1:end-1, :, :, :);
% kspace_info.DCF = kspace_info.DCF(2:end, :);
    
viewOrder = kspace_info.viewOrder;

if exist('recon_arms', 'var')
    if isnumeric(recon_arms)
        kspace = kspace(:, recon_arms, :, :);
        viewOrder = viewOrder(recon_arms);
    else
        kspace = kspace(:, 151:end, :, :);
        viewOrder = viewOrder(151:end);
    end
end

GA_steps = size(kx, 2);
Narms_total = size(kspace, 2);
Nframes = floor(Narms_total / Narms_per_frame);
Narms_total = Nframes * Narms_per_frame;
Ncoil = size(kspace, 4);
Nsample = size(kspace, 1);

kx = repmat(kx, [1, ceil(Narms_total / GA_steps)]);
ky = repmat(ky, [1, ceil(Narms_total / GA_steps)]);

kspace(:, Narms_total + 1 : end, :, :) = [];
viewOrder(Narms_total + 1 : end) = [];

kx = kx(:, viewOrder);
ky = ky(:, viewOrder);

kspace = reshape(kspace, [Nsample, Narms_per_frame, Nframes, Ncoil]);

Nsample_k = size(kx, 1);
kx = reshape(kx, [Nsample_k, Narms_per_frame, Nframes]);
ky = reshape(ky, [Nsample_k, Narms_per_frame, Nframes]);
    
%% DCF
% kspace_info.DCF(kspace_info.DCF > 0.2) = 0.2;
% kspace_info.DCF = kspace_info.DCF * 5;
    
%% spiral out
delay = 0;

kx_in = kx((1:floor(Nsample/2)) + delay, :, :);
ky_in = ky((1:floor(Nsample/2)) + delay, :, :);

DCF = kspace_info.DCF((1:floor(Nsample/2)) + delay, :);

Data.N = NUFFT.init(kx_in, ky_in, 1, [6, 6], matrix_size(1), matrix_size(1));
Data.N.W = DCF(:, 1);

Data.kSpace = kspace(1:floor(Nsample/2), :, :, :);
Data.first_est = NUFFT.NUFFT_adj(Data.kSpace, Data.N);
Data.sens_map = get_sens_map(Data.first_est, '2D');
Data.first_est = sum(Data.first_est .* conj(Data.sens_map), 4);

scale = max(abs(Data.first_est(:)));

para.Recon.no_comp = Ncoil;
para.Recon.weight_tTV = scale * para.weight_tTV; % temporal regularization weight
para.Recon.weight_sTV = scale * para.weight_sTV; % spatial regularization weight

[Image_2, para] = STCR_conjugate_gradient(Data, para);
%% save
save(sprintf('../../Dynamic_MaxGIRF/recon_data/%s_stcr_ttv_%03g_stv_%03g.mat', filename(1:end-4), para.weight_tTV, para.weight_sTV), 'Image_1', 'Image_2', 'para')
