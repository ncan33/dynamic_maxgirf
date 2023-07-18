
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

for i = 1 % loop through all dataset
    
    %% read data
    file_name = fullfile(all_dat(i).folder, all_dat(i).name);
    fprintf('--------------------------------- \n')
    fprintf('reading data:  %d/%d \n', i, nfile)
    fprintf('file name:     %s \n', file_name)
    fprintf('\n')
    
    load(file_name)
    
    narm    = size(kspace, 2);
    nview   = kspace_info.user_interleaves;
    
    %% reconstruction
    fprintf('gridding reconstruction \n')
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
    
    view_order_echo_1 = view_order(echo_idx == 1);
    view_order_echo_2 = view_order(echo_idx == 2);
    view_order_echo_1(narm_total + 1 : end) = [];
    view_order_echo_2(narm_total + 1 : end) = [];
    view_order_echo_1  = reshape(view_order_echo_1, [narm_frame, nframes]);
    view_order_echo_2  = reshape(view_order_echo_2, [narm_frame, nframes]);

    %% gridding 
    kx = kspace_info.kx_GIRF;
    ky = kspace_info.ky_GIRF;
    
    res = [kspace_info.user_ResolutionX, kspace_info.user_ResolutionY];
    fov = [kspace_info.user_FieldOfViewX, kspace_info.user_FieldOfViewY];
    matrix_size_keep    = round(fov ./ res / 2) * 2;
    matrix_size         = round(matrix_size_keep * 2 / 2) * 1.5;
    
    %% gridding recon
    N   = NUFFT.init(kx * matrix_size(1), ky * matrix_size(2), 1, [6, 6], matrix_size(1), matrix_size(1));
    N.W = kspace_info.DCF(:, 1);
    
    im_echo_1 = zeros([matrix_size, nframes, ncoil], 'single');
    im_echo_2 = zeros([matrix_size, nframes, ncoil], 'single');
    for iframe = 1:nframes
        
        idx_temp = view_order_echo_1(:, iframe);
        kspace_temp = zeros([nsample, nview, 1, ncoil]);
        kspace_temp(:, idx_temp, :, :) = kspace_echo_1(:, :, iframe, :);
        im_echo_1(:, :, iframe, :) = NUFFT.NUFFT_adj(kspace_temp, N);
        
        idx_temp = view_order_echo_2(:, iframe);
        idx_temp = idx_temp - nview;
        kspace_temp = zeros([nsample, nview, 1, ncoil]);
        kspace_temp(:, idx_temp, :, :) = kspace_echo_2(:, :, iframe, :);
        im_echo_2(:, :, iframe, :) = NUFFT.NUFFT_adj(kspace_temp, N);
        
    end
    im_echo_1 = squeeze(crop_half_FOV(im_echo_1, matrix_size_keep));
    im_echo_1 = fliplr(rot90(im_echo_1, -1));
    
    im_echo_2 = squeeze(crop_half_FOV(im_echo_2, matrix_size_keep));
    im_echo_2 = fliplr(rot90(im_echo_2, -1));
    
    %% save 
    f = imageMRI(sos([im_echo_1(:, :, end, :), im_echo_2(:, :, end, :)]));
    f = imageMRI(angle([im_echo_1(:, :, end, 1), im_echo_2(:, :, end, 1)]));
    if ifsave
        save_name   = sprintf('./recon_data/%s_recon.mat', all_dat(i).name(1:end-8));
%         figure_name = sprintf('./figure/%s_recon_naverage_%g.png', all_dat(i).name(1:end-8), naverage);
%         saveas(f, figure_name);
        save(save_name, 'im_echo_1', 'im_echo_2', 'kspace_info', '-v7.3');
    end

    toc
    
end
% close all