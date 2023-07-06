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
