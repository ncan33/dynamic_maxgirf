function recon_spiral_oioi(para)
   %% read data 
   
    load(para.raw_file)
    t1 = tic;
    if ~exist('kspace_info')
        kspace_info = kSpace_info;
        clear kSpace_info
    end
    if ~exist('kspace')
        kspace = kSpace;
        clear kSpace
    end
    para.kspace_info = kspace_info;
    
    TR = kspace_info.user_TR / 1000;
    para.Recon.narm = floor(para.Recon.temporal_resolution / TR);
    Narms_per_frame = para.Recon.narm;
    
    res = [kspace_info.user_ResolutionX, kspace_info.user_ResolutionY];
    
    kspace = permute(kspace, [1, 2, 4, 3]);
    
    matrix_size = round(para.Recon.FOV_recon ./ res / 2) * 2;
    
    para.Recon.image_size = matrix_size;
    
    matrix_size_keep = [kspace_info.user_FieldOfViewX, kspace_info.user_FieldOfViewX] ./ res;
    para.Recon.matrix_size_keep = round(matrix_size_keep);
    
    kx = kspace_info.kx_GIRF * matrix_size(1);
    ky = kspace_info.ky_GIRF * matrix_size(2);

    %%
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
    kspace_info.DCF(kspace_info.DCF > 0.2) = 0.2;
    kspace_info.DCF = kspace_info.DCF * 5;
    
    %% spiral out 1
    delay = 0;

    kx_in = kx((1:floor(Nsample/4)) + delay, :, :);
    ky_in = ky((1:floor(Nsample/4)) + delay, :, :);

    DCF = kspace_info.DCF((1:floor(Nsample/4)) + delay, :);
    
    Data.N = NUFFT.init(kx_in, ky_in, 1, [6, 6], matrix_size(1), matrix_size(1));
    Data.N.W = DCF(:, 1);
    
    Data.kSpace = kspace(1:floor(Nsample/4), :, :, :);
    Data.first_est = NUFFT.NUFFT_adj(Data.kSpace, Data.N);
    Data.sens_map = get_sens_map(Data.first_est, '2D');
    Data.first_est = sum(Data.first_est .* conj(Data.sens_map), 4);
    
    scale = max(abs(Data.first_est(:)));
    
    para.Recon.no_comp = Ncoil;
    para.Recon.weight_tTV = scale * para.weight_tTV; % temporal regularization weight
    para.Recon.weight_sTV = scale * para.weight_sTV; % spatial regularization weight

    [image_out_1, para] = STCR_conjugate_gradient(Data, para);
    para.cost_spiral_out_1 = para.Cost;
    para.cpu_time_spiral_out_1 = para.CPUtime;
    para = rmfield(para, 'Cost');
    para = rmfield(para, 'CPUtime');
    para.Recon.step_size = 2;
    
    %% spiral in 1
    delay = 0;
    
    kx_out = kx(floor(Nsample/4)+1:Nsample/2 + delay, :, :);
    ky_out = ky(floor(Nsample/4)+1:Nsample/2 + delay, :, :);

    DCF = kspace_info.DCF(floor(Nsample/4)+1:Nsample/2 + delay, :);
    
    Data.N = NUFFT.init(kx_out, ky_out, 1, [6, 6], matrix_size(1), matrix_size(1));
    Data.N.W = DCF(:, 1);
    
    Data.kSpace = kspace(floor(Nsample/4)+1:Nsample/2 + delay, :, :, :);
    Data.first_est = NUFFT.NUFFT_adj(Data.kSpace, Data.N);
    Data.first_est = sum(Data.first_est .* conj(Data.sens_map), 4);

    [image_in_1, para] = STCR_conjugate_gradient(Data, para);
    para.cost_spiral_in_1 = para.Cost;
    para.cpu_time_spiral_in_1 = para.CPUtime;
    para = rmfield(para, 'Cost');
    para = rmfield(para, 'CPUtime');
    para.Recon.step_size = 2;
    
    %% spiral out 2
    delay = 0;

    kx_in = kx((Nsample/2 + 1 : Nsample/2 + floor(Nsample/4)) + delay, :, :);
    ky_in = ky((Nsample/2 + 1 : Nsample/2 + floor(Nsample/4)) + delay, :, :);

    DCF = kspace_info.DCF((Nsample/2 + 1 : Nsample/2 + floor(Nsample/4)) + delay, :);
    
    Data.N = NUFFT.init(kx_in, ky_in, 1, [6, 6], matrix_size(1), matrix_size(1));
    Data.N.W = DCF(:, 1);
    
    Data.kSpace = kspace((Nsample/2 + 1 : Nsample/2 + floor(Nsample/4)) + delay, :, :, :);
    Data.first_est = NUFFT.NUFFT_adj(Data.kSpace, Data.N);
    Data.first_est = sum(Data.first_est .* conj(Data.sens_map), 4);

    [image_out_2, para] = STCR_conjugate_gradient(Data, para);
    para.cost_spiral_out_2 = para.Cost;
    para.cpu_time_spiral_out_2 = para.CPUtime;
    para = rmfield(para, 'Cost');
    para = rmfield(para, 'CPUtime');
    para.Recon.step_size = 2;
    
    %% spiral in 2
    delay = 0;

    kx_in = kx(Nsample/2 + floor(Nsample/4) + 1 + delay : Nsample, :, :);
    ky_in = ky(Nsample/2 + floor(Nsample/4) + 1 + delay : Nsample, :, :);

    DCF = kspace_info.DCF(Nsample/2 + floor(Nsample/4) + 1 + delay : Nsample, :);
    
    Data.N = NUFFT.init(kx_in, ky_in, 1, [6, 6], matrix_size(1), matrix_size(1));
    Data.N.W = DCF(:, 1);
    
    Data.kSpace = kspace(Nsample/2 + floor(Nsample/4) + 1 + delay : Nsample, :, :, :);
    Data.first_est = NUFFT.NUFFT_adj(Data.kSpace, Data.N);
    Data.first_est = sum(Data.first_est .* conj(Data.sens_map), 4);

    [image_in_2, para] = STCR_conjugate_gradient(Data, para);
    para.cost_spiral_in_2 = para.Cost;
    para.cpu_time_spiral_in_2 = para.CPUtime;
    para = rmfield(para, 'Cost');
    para = rmfield(para, 'CPUtime');
    para.Recon.step_size = 2;
    
    %% save
    if ~isfolder([para.vol_num, '/recon_data'])
        mkdir([para.vol_num, '/recon_data'])
    end
    para.recon_time = toc(t1);
    save([para.vol_num, '/recon_data/', para.file_name(1:end-4), '_stcr.mat'], 'image_in_1', 'image_out_1', 'image_in_2', 'image_out_2', 'para')

end
