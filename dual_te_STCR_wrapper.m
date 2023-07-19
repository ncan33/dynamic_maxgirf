function [im_echo_1, im_echo_2, NUFFT_im_echo_1, NUFFT_im_echo_2, ...
    kspace_info, para] = dual_te_STCR_wrapper(narm_frame, tTV, sTV, ...
    niter, ifsave, ifGPU, ifplot, path)
    % This is a wrapper function for spatiotemporally constrained reconstruction
    % on dual TE variable density spiral raw RTHawk data. It calls STCR on my
    % the data of this repository
    
    arguments
        narm_frame
        tTV
        sTV
        niter
        ifsave = 0
        ifGPU = 1
        ifplot = 1;
        path = '/server/sdata/ncan/mri_data/disc/lung/vol0457_20221021/raw_hawk/usc_disc_yt_2022_10_21_133643_dual-te_dynamic.mat'
    end
    
    %% add ismrmd and mfile path
    addpath ./util/mfile/functions/
    addpath ./util/mfile/registrtation/
    addpath ./util/mfile/quantification/
    addpath ./util/mfile/vdspiral/
    addpath ./util/

    ccc %check that setup.m is run
    if ~isfolder('./recon_data')
        mkdir ./recon_data
    end
    
    all_dat = dir(path);
    %all_dat     = dir('/server/sdata/ncan/mri_data/disc/lung/vol0457_20221021/raw_hawk/usc_disc_yt_2022_10_21_133643_dual-te_dynamic.mat');
    nfile       = length(all_dat);
    % naverage    = 1;
    % nprep       = 5;
    file_index = 1;

    %% read data
    file_name = fullfile(all_dat(file_index).folder, all_dat(file_index).name);
    fprintf('--------------------------------- \n')
    fprintf('reading data:  %d/%d \n', file_index, nfile)
    fprintf('file name:     %s \n', file_name)
    fprintf('\n')

    load(file_name)

    narm = size(kspace, 2); % total number of spiral arms acquired during scan
    nview = kspace_info.user_interleaves; % total number of interleaves in original kspace trajectory

    %% get echo index
    tic

    view_order = kspace_info.viewOrder; % view order

    echo_idx = view_order > nview;
    echo_idx = echo_idx + 1; % if the element in echo_idx == 1, then echo_1,
                             % and when echo_idx == 2, then echo_2

    %% size
    nsample = size(kspace, 1); % number of samples
    ncoil = size(kspace, 3); % number of coils

    %% set reconstruction parameters
    para.Recon.narm = narm_frame; % number of spiral arms per frame
    para.Recon.time_frames = 'all'; % set to 'all' for reconstructructing all frames. specify certain range if you wish, i.e., 1:100

    para.weight_tTV = tTV; % expected range is 1e-1 to 1e-5 
    para.weight_sTV = sTV; % expected range is 1e-1 to 1e-6

    para.setting.ifplot = ifplot; % display image and cost during reconstruction
    para.setting.ifGPU = ifGPU; % set to 1 when you want to use GPU

    para.Recon.no_comp = ncoil;     % number of coils
    %para.Recon.FOV = 1.25;         % reconstruction FOV

    para.Recon.epsilon = eps('single');
    para.Recon.step_size = 2;
    para.Recon.ifContinue = 0;
    para.Recon.noi = niter; % number of iterations
    para.Recon.type = '2D Spiral server';
    para.Recon.break = 1;

    %% split kspace
    kspace_echo_1 = kspace(:, echo_idx == 1, :); % raw kspace data for echo 1
    kspace_echo_2 = kspace(:, echo_idx == 2, :); % raw kspace data for echo 2
    %clear kspace

    narm_total = min(size(kspace_echo_1, 2), size(kspace_echo_2, 2)); % narm after splitting kspace

    %% orgnize the data to frames
    nframes = floor(narm_total / narm_frame); % I don't know if this is necessary
    narm_total = nframes * narm_frame; % I don't know if this is necessary

    kspace_echo_1(:, narm_total + 1 : end, :) = []; % discard excess kspace
    kspace_echo_2(:, narm_total + 1 : end, :) = []; % discard excess kspace

    kspace_echo_1  = reshape(kspace_echo_1, [nsample, narm_frame, nframes, ncoil]);
    kspace_echo_2  = reshape(kspace_echo_2, [nsample, narm_frame, nframes, ncoil]);

    %% kspace trajectory and view order
    view_order_echo_1 = view_order(echo_idx == 1); % view_order for echo_1
    view_order_echo_2 = view_order(echo_idx == 2); % view_order for echo_2
    view_order_echo_1(narm_total + 1 : end) = []; % discard excess views
    view_order_echo_2(narm_total + 1 : end) = []; % discard excess views
    view_order_echo_1  = reshape(view_order_echo_1, [narm_frame, nframes]);
    view_order_echo_2  = reshape(view_order_echo_2, [narm_frame, nframes]);

    kx = kspace_info.kx_GIRF; % kx
    ky = kspace_info.ky_GIRF; % ky

    kx_echo_1 = zeros(nsample, narm_frame, nframes); %kx_echo_1 = zeros(size(kspace_echo_1,1), size(kspace_echo_1,2), size(kspace_echo_1,3));
    kx_echo_2 = kx_echo_1;
    ky_echo_1 = kx_echo_1;
    ky_echo_2 = kx_echo_1;

    for i = 1:narm_frame
        for j = 1:nframes
            kx_echo_1(:,i,j) = kx(:, view_order_echo_1(i,j));
            kx_echo_2(:,i,j) = kx(:, view_order_echo_2(i,j) - 10);

            ky_echo_1(:,i,j) = ky(:, view_order_echo_1(i,j));
            ky_echo_2(:,i,j) = ky(:, view_order_echo_2(i,j) - 10);
        end
    end

    [im_echo_1, NUFFT_im_echo_1, para] = STCR(kspace_info, kspace_echo_1, ...
        kx_echo_1, ky_echo_1, para);

    %% perform STCR for echo 2
    [im_echo_2, NUFFT_im_echo_2, para] = STCR(kspace_info, kspace_echo_2, ...
        kx_echo_2, ky_echo_2, para);

    %% save 
    f_mag = imageMRI(abs([im_echo_1(:,:,end), im_echo_2(:,:,end)]));
    %f_phase = imageMRI(angle(im_echo_1(:,:,end), im_echo_2(:,:,end)));
    if ifsave
        save_name   = sprintf(['./recon_data/',num2str(narm_frame),'arm_',num2str(para.weight_tTV),'_tTV_',num2str(para.weight_sTV),'_sTV_','%s_recon.mat'], all_dat(file_index).name(1:end-8));
    %         figure_name = sprintf('./figure/%s_recon_naverage_%g.png', all_dat(i).name(1:end-8), naverage);
    %         saveas(f, figure_name);
        save(save_name, 'im_echo_1', 'im_echo_2', 'NUFFT_im_echo_1', 'NUFFT_im_echo_2', 'kspace_info', '-v7.3');
    end

    toc
end
