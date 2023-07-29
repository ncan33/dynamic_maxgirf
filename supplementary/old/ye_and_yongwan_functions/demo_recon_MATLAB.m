clear; close all; clc;

%% add ismrmd and mfile path
%addpath /server/home/ncan/ismrmrd % hermia
addpath /Users/nick/research/mrel/ismrmrd/matlab % mac
addpath ./matlab/

%% load raw data and save recon file directories
para.parent_dir = '/Volumes/SamsungSSD/Research/mrel/mri_data/speech';
para.dir.subj = 'sub001';
para.dir.task = '2drt_01_vcv1_r1';
para.dir.raw_file = fullfile(para.parent_dir, para.dir.subj, '2drt/raw', sprintf('%s_%s_raw.h5', para.dir.subj, para.dir.task));
para.dir.save_recon = fullfile(para.parent_dir, para.dir.subj, 'recon', sprintf('%s_%s_recon.h5', para.dir.subj, para.dir.task));
if ~isfolder(fullfile(para.parent_dir, para.dir.subj, 'recon'))
    mkdir(fullfile(para.parent_dir, para.dir.subj, 'recon'))
end

%% set parameters
para.setting.ifplot = 1;        % display image and cost during reconstruction
para.setting.ifGPU = 0;         % set to 1 when you want to use GPU
para.Recon.time_frames = 1:100; % set to 'all' for reconstructing all frames

para.weight_tTV = 0.08;         % temporal regularization parameter
para.weight_sTV = 0;            % spatial regularization parameter

para.Recon.narm = 2;            % number of arms per frame
para.Recon.FOV = 1.25;          % reconstruction FOV

para.Recon.epsilon = eps('single');
para.Recon.step_size = 2;
para.Recon.ifContinue = 0;
para.Recon.noi = 150; % number of iterations
para.Recon.type = '2D Spiral server';
para.Recon.break = 1;

%% do the recon
%speech_open_dataset_reconstruction(para)

%% read data
dset = ismrmrd.Dataset(para.dir.raw_file, 'dataset');
data = dset.readAcquisition();
head = ismrmrd.xml.deserialize(dset.readxml);

matrix_size = [head.encoding.reconSpace.matrixSize.x, head.encoding.reconSpace.matrixSize.y];
para.Recon.matrix_size = matrix_size;

%% fix dimensions for kspace, kx and ky
kSpace = cat(3, data.data{:});
kSpace = permute(kSpace, [1, 3, 2]);
kSpace = kSpace / max(abs(kSpace(:))) * 1e6; % scale k-space
disp(['size(kSpace):', newline, num2str(size(kSpace))])

ktraj = cat(3, data.traj{:});
kx = squeeze(ktraj(1, :, :));
ky = squeeze(ktraj(2, :, :));

kx = kx * matrix_size(1);
ky = ky * matrix_size(2);
disp(['size(kx):', newline, num2str(size(kx))])
disp(['size(ky):', newline, num2str(size(ky))])

w = squeeze(ktraj(3, :, 1)).'; % density compensation 
disp(['size(w):', newline, num2str(size(ky))])

%clearvars -except kSpace kx ky w para

[sx, ns, nc] = size(kSpace);

narm = para.Recon.narm;

nof = floor(ns/narm);
kSpace(:, nof*narm+1:end, :) = [];
kx(:, nof*narm+1:end) = [];
ky(:, nof*narm+1:end) = [];
disp(['size(kx):', newline, num2str(size(kx))])
disp(['size(ky):', newline, num2str(size(ky))])

kSpace = reshape(kSpace, [sx, narm, nof, nc]);
kx = reshape(kx, [sx, narm, nof]);
ky = reshape(ky, [sx, narm, nof]);
disp(['size(kSpace):', newline, num2str(size(kSpace))])
disp(['size(kx):', newline, num2str(size(kx))])
disp(['size(ky):', newline, num2str(size(ky))])

%% select time frames
if isnumeric(para.Recon.time_frames)
    time_frames = para.Recon.time_frames;
    
    kx = kx(:, :, time_frames);
    ky = ky(:, :, time_frames);
    kSpace = kSpace(:, :, time_frames, :);
end

Data.N = NUFFT.init(kx*para.Recon.FOV, ky*para.Recon.FOV, 1, [4, 4], para.Recon.matrix_size(1)*para.Recon.FOV, para.Recon.matrix_size(1)*para.Recon.FOV);
Data.N.W = w;

Data.kSpace = kSpace;
Data.first_est = NUFFT.NUFFT_adj(Data.kSpace, Data.N);

scale = max(abs(Data.first_est(:)));

Data.sens_map = get_sens_map(Data.first_est, '2D');
Data.first_est = sum(Data.first_est .* conj(Data.sens_map), 4);

%% set parameters
para.Recon.no_comp = nc;
para.Recon.weight_tTV = scale * para.weight_tTV; % temporal regularization weight
para.Recon.weight_sTV = scale * para.weight_sTV; % spatial regularization weight

%clearvars -except Data para


%% conjugate gradient reconstruction
[Image_recon, para] = STCR_conjugate_gradient(Data, para);
Image_recon = rot90(abs(Image_recon));
Image_recon = crop_half_FOV(Image_recon, para.Recon.matrix_size);
%{
%% save reconstruction
if isfile(para.dir.save_recon)
    delete(para.dir.save_recon)
end
h5create(para.dir.save_recon,'/recon', size(Image_recon))
h5write(para.dir.save_recon, '/recon', Image_recon)
h5create(para.dir.save_recon,'/lambda_t', [1, 1])
h5write(para.dir.save_recon, '/lambda_t', para.weight_tTV)
h5create(para.dir.save_recon,'/lambda_s', [1, 1])
h5write(para.dir.save_recon, '/lambda_s', para.weight_sTV)
h5create(para.dir.save_recon,'/iterations', [1, 1])
h5write(para.dir.save_recon, '/iterations', para.Recon.noi)
end
%}