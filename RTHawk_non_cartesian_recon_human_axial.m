% demo_non_cartesian_recon_human_axial.m
% Written by Nejat Can
% Email: ncan@usc.edu
% Started: 08/21/2023

%% Clean slate
close all; clear all; clc;

%% Set source directories
src_directory = '/server/home/ncan/GitHub/lowfield_maxgirf';
ismrmrd_directory = '/server/home/ncan/ismrmrd';

%% Add source directories to search path
rmpath('./util/functions')
rmpath('./util/mfile/vdspiral')
rmpath('./util/mfile/quantification')
rmpath('./util/mfile/registrtation')
rmpath('./util/mfile/functions')

addpath(genpath(src_directory));

%% Define data directory
data_path = '/server/sdata/ncan/mri_data/disc/lung/vol0457_20221021/raw_hawk/usc_disc_yt_2022_10_21_133643_dual-te_dynamic.mat';
B0map_fullpath = '/server/home/ncan/GitHub/dynamic_maxgirf/B0map_zeros_420_420_59.mat';

%% Set reconstruction options
% "phase_sign" and "read_sign" can be determined only from Siemens raw data 
% format now until the ISMRMRD format includes these as part of its header
user_opts.narm_frame           = 5;                % number of arms per frame
user_opts.vds_factor           = 75;
user_opts.discard_pre          = 20;
user_opts.discard_post         = 20;
user_opts.N1                   = 420;              % reconstruction matrix size along the row direction
user_opts.N2                   = 420;              % reconstruction matrix size along the column direction
user_opts.max_iterations       = 45;               % maximum number of LSQR iterations
user_opts.tol                  = 1e-5;             % LSQR tolerance
user_opts.static_B0_correction = 1;                % static off-resonance correction: 1=yes, 0=no
user_opts.Lmax                 = 20;               % maximum rank of the SVD approximation of a higher-order encoding matrix
user_opts.L                    = 5;                % rank of the SVD approximation of a higher-order encoding matrix
user_opts.support_constraint   = 1;                % sse a support constraint using a voxel mask created with ESPIRiT

%% Define an output filename
[filepath,filename,ext] = fileparts(data_path);
output_filename = sprintf('%s_%dx%d', filename, user_opts.N1, user_opts.N2);

%% Load a static off-resonance map [Hz]
load(B0map_fullpath);
B0map = B0map_nlinv;

if 0
%% Perform NUFFT reconstruction
[im_nufft, header_nufft, r_dcs_nufft] = siemens_gridding_recon(ismrmrd_noise_path, ismrmrd_data_path, siemens_dat_path, user_opts);
save(sprintf('%s_nufft', output_filename), 'im_nufft', 'header_nufft', 'r_dcs_nufft', 'user_opts', '-v7.3');

%% Perform SENSE reconstruction
[im_sense, header_sense, r_dcs_sense, output_sense] = siemens_sense_recon(ismrmrd_noise_path, ismrmrd_data_path, siemens_dat_path, user_opts);
save(sprintf('%s_sense', output_filename), 'im_sense', 'header_sense', 'r_dcs_sense', 'output_sense', 'user_opts', '-v7.3');

%% Perform CP-based MaxGIRF reconstruction
[im_cpr, header_cpr, r_dcs_cpr, output_cpr] = siemens_maxgirf_cp_recon(ismrmrd_noise_path, ismrmrd_data_path, siemens_dat_path, B0map_nlinv, user_opts);
save(sprintf('%s_cpr', output_filename), 'im_cpr', 'header_cpr', 'r_dcs_cpr', 'output_cpr', 'user_opts', '-v7.3');

%% Perform CG-based MaxGIRF reconstruction
[im_maxgirf, header_maxgirf, r_dcs_maxgirf, output_maxgirf] = siemens_maxgirf_cg_recon(ismrmrd_noise_path, ismrmrd_data_path, siemens_dat_path, B0map_nlinv, user_opts);
save(sprintf('%s_maxgirf', output_filename), 'im_maxgirf', 'header_maxgirf', 'r_dcs_maxgirf', 'output_maxgirf', 'user_opts', '-v7.3');
end

%% Perform NUFFT reconstruction (single-GPU)
[im_nufft_gpu, header_nufft, r_dcs_nufft] = RTHawk_gridding_recon_gpu(data_path, user_opts);
save(sprintf('%s_nufft_gpu', output_filename), 'im_nufft_gpu', 'header_nufft', 'r_dcs_nufft', 'user_opts', '-v7.3');

if 0
%% Perform CG-based MaxGIRF reconstruction (single-GPU)
[im_maxgirf_gpu, header_maxgirf, r_dcs_maxgirf, output_maxgirf] = siemens_maxgirf_cg_recon_single_gpu(ismrmrd_noise_path, ismrmrd_data_path, siemens_dat_path, B0map_nlinv, user_opts);
save(sprintf('%s_maxgirf_single_gpu_supp%d_iter%d', output_filename, user_opts.support_constraint, user_opts.max_iterations), 'im_maxgirf_gpu', 'header_maxgirf', 'r_dcs_maxgirf', 'output_maxgirf', 'user_opts', '-v7.3');

%% Perform King's method reconstruction
[im_king, header_king, r_dcs_king] = siemens_king_method_recon(ismrmrd_noise_path, ismrmrd_data_path, siemens_dat_path, user_opts);
save(sprintf('%s_king', output_filename), 'im_king', 'header_king', 'r_dcs_king', 'user_opts', '-v7.3');

%% Perform SENSE reconstruction (single-GPU)
[im_sense_gpu, header_sense, r_dcs_sense, output_sense] = siemens_sense_recon_gpu(ismrmrd_noise_path, ismrmrd_data_path, siemens_dat_path, user_opts);
save(sprintf('%s_sense_gpu', output_filename), 'im_sense_gpu', 'header_sense', 'r_dcs_sense', 'output_sense', 'user_opts', '-v7.3');
end

if 0
%% Perform CP-based MaxGIRF reconstruction (single-GPU)
[im_cpr_gpu, header_cpr, r_dcs_cpr, output_cpr] = siemens_maxgirf_cp_recon_gpu(ismrmrd_noise_path, ismrmrd_data_path, siemens_dat_path, B0map_nlinv, user_opts);
save(sprintf('%s_cpr_gpu', output_filename), 'im_cpr_gpu', 'header_cpr', 'r_dcs_cpr', 'output_cpr', 'user_opts', '-v7.3');

%% Perform CG-based MaxGIRF reconstruction (single-GPU)
[im_maxgirf_gpu, header_maxgirf, r_dcs_maxgirf, output_maxgirf] = siemens_maxgirf_cg_recon_single_gpu(ismrmrd_noise_path, ismrmrd_data_path, siemens_dat_path, B0map_nlinv, user_opts);
save(sprintf('%s_maxgirf_multi_gpu_supp%d_iter%d', output_filename, user_opts.support_constraint, user_opts.max_iterations), 'im_maxgirf_gpu', 'header_maxgirf', 'r_dcs_maxgirf', 'output_maxgirf', 'user_opts', '-v7.3');

%% Perform CG-based MaxGIRF reconstruction (multi-GPU)
[im_maxgirf_gpu, header_maxgirf, r_dcs_maxgirf, output_maxgirf] = siemens_maxgirf_cg_recon_multi_gpu(ismrmrd_noise_path, ismrmrd_data_path, siemens_dat_path, B0map_nlinv, user_opts);
save(sprintf('%s_maxgirf_multi_gpu_supp%d_iter%d', output_filename, user_opts.support_constraint, user_opts.max_iterations), 'im_maxgirf_gpu', 'header_maxgirf', 'r_dcs_maxgirf', 'output_maxgirf', 'user_opts', '-v7.3');
end
