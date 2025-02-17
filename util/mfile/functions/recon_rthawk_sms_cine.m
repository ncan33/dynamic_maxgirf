function recon_rthawk_sms_cine(vol_number)

area = 'heart';
n_vol = length(vol_number);

for ivol = 1:n_vol
    path_server = dir(sprintf('/mnt/sdata_new/mri_data/disc/%s/vol%04g_*', area, vol_number(ivol)));
    path_server = fullfile(path_server.folder, path_server.name);
    
    path_yt = strfind(path_server, 'mri_data');
    path_yt = path_server(path_yt:end);
    path_yt = ['/mnt/sdata_new/ytian/', path_yt];
    
    path_yt_raw = [path_yt, '/raw_hawk/'];
    if ~isfolder(path_yt_raw)
        mkdir(path_yt_raw)
    end
    
    path_yt_recon = [path_yt, '/recon_data/new_6/'];
    if ~isfolder(path_yt_recon)
        mkdir(path_yt_recon)
    end
    
    all_mat = dir([path_yt_raw, 'usc_disc_yt*sms*.mat']);
    nfile = length(all_mat);
    
    %%
    FOV_recon = [500, 500]; % mm
    
    %% set settings
    para.setting.ifplot         = 0;        % plot convergence during reconstruction
    para.setting.ifGPU          = 1;        % set to 1 when you want to use GPU
    
    %% set recon parameters
    para.weight_tTV = 0.002;                      % temporal TV regularizaiton parameter (normalized by F^T d)
    para.weight_sTV = 0.0005;                        % spatial TV regularizaiton parameter (normalized by F^T d)
    
    para.Recon.epsilon      = eps('single');        % small vale to avoid singularity in TV constraint
    para.Recon.step_size    = 2;                    % initial step size
    para.Recon.noi          = 40;                  % number of CG iterations
    para.Recon.type         = '2D spiral sms server'; % stack of spiral
    para.Recon.break        = 1;                    % stop iteration if creteria met. Otherwise will run to noi
    
    for ifile = 2:nfile
        %% read data
        fprintf('reconstructing data: %s vol%04g %s \n',area, vol_number(ivol), all_mat(ifile).name)
        file_name = fullfile(all_mat(ifile).folder, all_mat(ifile).name);
        load(file_name)
        
        n_tr_per_frame_per_rr = 4;
        
        minimal_n_tr_per_frame = n_tr_per_frame_per_rr * max(kSpace_info.ViewSegmentIndex + 1);
        
        %% load parameters
        nsample = kSpace_info.extent(1);
        ncoil   = kSpace_info.extent(2);
        nsms    = max(kSpace_info.RFIndex) + 1;
        res     = [kSpace_info.user_ResolutionX, kSpace_info.user_ResolutionY];
        
        matrix_size = round(FOV_recon ./ res / 2) * 2;
        para.Recon.image_size = matrix_size;
        matrix_size_keep = [kSpace_info.user_FieldOfViewX, kSpace_info.user_FieldOfViewX] ./ res;
        para.Recon.matrix_size_keep = round(matrix_size_keep);
        para.Recon.no_comp = ncoil;
        
        n_tr_per_frame_per_rr_ = n_tr_per_frame_per_rr;
        minimal_n_tr_per_frame_ = minimal_n_tr_per_frame;
        
        %% find trigger position
        if isfield(kSpace_info, 'ViewSegmentIndex')
            trig_pos = find(diff(kSpace_info.ViewSegmentIndex) == 1) + 1;
            %         trig_pos(end) = [];
            trig_pos = [1, trig_pos];
            %         trig_pos = [trig_pos, size(kSpace, 2)];
        else
            trig_pos = local_max(kSpace_info.TimeSinceTrig);
            trig_pos = trig_pos + 1;
        end
        
        if trig_pos(end) > length(kSpace_info.TimeSinceTrig)
            trig_pos(end) = [];
        end
        
        if trig_pos(end) ~= length(kSpace_info.TimeSinceTrig)
            trig_pos(end + 1) = length(kSpace_info.TimeSinceTrig);
        end
        
        trig_idx = [];
        for i = 1:length(trig_pos) - 1
            n_tr_since_trig_temp = 1 : trig_pos(i + 1) - trig_pos(i);
            trig_idx = [trig_idx, n_tr_since_trig_temp];
        end
        trig_idx(end+1) = trig_idx(end) + 1;
        
        cardiac_phase_idx = ceil(trig_idx / n_tr_per_frame_per_rr_);
        
        %% YT 2021/11/17
        n_cardiac_phase = floor(median(cardiac_phase_idx([trig_pos(2:end-1) - 1, trig_pos(end)])));
        n_tr_per_phase = sum(cardiac_phase_idx == 1);
        
        %% remove data
        view_order = kSpace_info.viewOrder;
        rf_idx = kSpace_info.RFIndex;
        idx_remove = cardiac_phase_idx > n_cardiac_phase;
        
        kSpace(:, idx_remove, :) = [];
        view_order(idx_remove) = [];
        rf_idx(idx_remove) = [];
        cardiac_phase_idx(idx_remove) = [];
        trig_idx(idx_remove) = [];
        
        trig_pos = find(trig_idx == 1);
        n_rr = length(trig_pos);
        trig_pos(end + 1) = length(view_order) + 1;
        
        n_tr_per_rr = n_cardiac_phase * n_tr_per_frame_per_rr;
        for irr = 1 : n_rr
            n_tr_this_rr = trig_pos(irr + 1) - trig_pos(irr);
            if n_tr_this_rr < n_tr_per_rr
                %             keyboard
                n_missing_tr = n_tr_per_rr - n_tr_this_rr;
                rr_idx = trig_pos(irr):(trig_pos(irr + 1)-1);
                idx_end = floor(length(rr_idx) / n_tr_per_frame_per_rr) * n_tr_per_frame_per_rr;
                
                idx_add_on = rr_idx(idx_end - n_tr_per_frame_per_rr + 1 : idx_end);
                idx_add_on = repmat(idx_add_on, [1, ceil(n_missing_tr / n_tr_per_frame_per_rr)]);
                idx_add_on = idx_add_on(end - n_missing_tr + 1 : end);
                %             idx_add_on = rr_idx(idx_end - n_missing_tr + 1 : idx_end);
                
                kSpace = cat(2, kSpace(:, 1:rr_idx(end), :), kSpace(:, idx_add_on, :), kSpace(:, rr_idx(end)+1:end, :));
                %             cardiac_phase_idx = cat(2, cardiac_phase_idx(1:rr_idx(end)),  , cardiac_phase_idx(rr_idx(end)+1:end));
                view_order = cat(2, view_order(1:rr_idx(end)), view_order(idx_add_on), view_order(rr_idx(end)+1:end));
                rf_idx = cat(2, rf_idx(1:rr_idx(end)), rf_idx(idx_add_on), rf_idx(rr_idx(end)+1:end));
                trig_pos(irr+1:end) = trig_pos(irr+1:end) + n_missing_tr;
                
            end
        end
        cardiac_phase_idx = vec(repmat(1:n_cardiac_phase, [n_tr_per_frame_per_rr, n_rr]))';
        
        %% sort data
        kspace_sort = zeros([nsample, n_tr_per_phase, n_cardiac_phase, ncoil], 'single');
        view_sort = ones([n_tr_per_phase, n_cardiac_phase], 'single');
        rf_sort = - ones([n_tr_per_phase, n_cardiac_phase], 'single');
        
        for iphase = 1:n_cardiac_phase
            idx = cardiac_phase_idx == iphase;
            n_tr_this_phase = sum(idx);
            idx_ = 1:n_tr_this_phase;
            
            kspace_sort(:, idx_, iphase, :) = kSpace(:, idx, :);
            view_sort(idx_, iphase) = view_order(idx);
            rf_sort(idx_, iphase) = rf_idx(idx);
        end
        
        %{
    n_cardiac_phase = max(cardiac_phase_idx);
    while sum(n_cardiac_phase == cardiac_phase_idx) < minimal_n_tr_per_frame_
        n_cardiac_phase = n_cardiac_phase - 1;
        if n_cardiac_phase == 0
            disp('wrong')
            break
        end
    end
    
    n_tr_per_phase = sum(cardiac_phase_idx == 1);
    
    %% remove first few k-space that are not in steady-state
%     kSpace_info.RFIndex(1:9) = -1;

    %% sort data
    kspace_sort = zeros([nsample, n_tr_per_phase, n_cardiac_phase, ncoil], 'single');
    view_sort = ones([n_tr_per_phase, n_cardiac_phase], 'single');
    rf_sort = - ones([n_tr_per_phase, n_cardiac_phase], 'single');
    
    for iphase = 1:n_cardiac_phase
        idx = cardiac_phase_idx == iphase;
        n_tr_this_phase = sum(idx);
        idx_ = 1:n_tr_this_phase;
        
        kspace_sort(:, idx_, iphase, :) = kSpace(:, idx, :);
        view_sort(idx_, iphase) = kSpace_info.viewOrder(idx);
        rf_sort(idx_, iphase) = kSpace_info.RFIndex(idx);
    end
    
%     kspace_sort(:, :, 1, :) = [];
%     view_sort(:, 1) = [];
%     rf_sort(:, 1) = [];
%     n_cardiac_phase = n_cardiac_phase - 1;
        %}
        nan_idx = rf_sort == -1;
        
        %% remove first 9 spirals that are not in steady-state, borrow data from a close cardaic phase
        
        n_cardiac_phase_cycle_1 = max(cardiac_phase_idx(kSpace_info.ViewSegmentIndex==0));
        if sum(cardiac_phase_idx(kSpace_info.ViewSegmentIndex==0) == n_cardiac_phase_cycle_1) < 6
            n_cardiac_phase_cycle_1 = n_cardiac_phase_cycle_1 - 1;
        end
        idx = find(cardiac_phase_idx == n_cardiac_phase_cycle_1);
        kspace_sort(:, 1:6, 1, :) = kSpace(:, idx(1:6), :);
        rf_sort(1:6, 1) = kSpace_info.RFIndex(idx(1:6));
        view_sort(1:6, 1) = kSpace_info.viewOrder(idx(1:6));
        kspace_sort(:, 1:3, 2, :) = kSpace(:, idx(1:3), :);
        rf_sort(1:3, 2) = kSpace_info.RFIndex(idx(1:3));
        view_sort(1:3, 2) = kSpace_info.viewOrder(idx(1:3));
        %{
        %}
        %% trajectory, phase, gridding
        kx = kSpace_info.kx_GIRF(:, view_sort) * matrix_size(1);
        ky = kSpace_info.ky_GIRF(:, view_sort) * matrix_size(2);
        
        kx = reshape(kx, [nsample, n_tr_per_phase, n_cardiac_phase]);
        ky = reshape(ky, [nsample, n_tr_per_phase, n_cardiac_phase]);
        
        Data.N = NUFFT.init(kx, ky, 1, [6, 6], matrix_size(1), matrix_size(1));
        Data.N.W = kSpace_info.DCF(:, 1);
        
        phase_index = rf_sort - 1;
        
        phase_mod = ([1:nsms] - ceil(nsms/2)) * 2 * pi / nsms;
        phase_mod = phase_index .* reshape(phase_mod, [1, 1, nsms]);
        phase_mod = exp(1i * phase_mod);
        phase_mod = phase_mod .* ~ nan_idx;
        phase_mod = permute(phase_mod, [4, 1, 2, 5, 3]);
        Data.phase_mod = phase_mod;
        
        %% preparation for iterative recon
        Data.kSpace = kspace_sort;
        Data.first_est = NUFFT.NUFFT_adj(kspace_sort .* conj(phase_mod), Data.N);
        Data.sens_map = get_sens_map(Data.first_est, 'SMS');
        Data.first_est = sum(Data.first_est .* conj(Data.sens_map), 4);
        Data.phase_mod = phase_mod;
        
        scale = max(abs(Data.first_est(:)));
        
        para.Recon.weight_tTV = scale * para.weight_tTV; % temporal regularization weight
        para.Recon.weight_sTV = scale * para.weight_sTV; % spatial regularization weight
        
        %% iterative stcr
        if strfind(file_name, 'sms_1')
            Image_recon = Data.first_est;
        else
            [Image_recon, para] = STCR_conjugate_gradient_circ(Data, para);
        end
        
        %% crop, rotate, save image
        Image_recon = fliplr(rot90(Image_recon, -1));
        Image_recon = abs(Image_recon);
        Image_recon = crop_half_FOV(Image_recon, para.Recon.matrix_size_keep);
        Image_recon = squeeze(Image_recon);
        
        if strfind(file_name, 'sms_1')
            save(sprintf('%s%s_gridding_frametime_%.2g.mat', path_yt_recon, all_mat(ifile).name(1:end-4), kSpace_info.user_TR * n_tr_per_frame_per_rr / 1000), 'Image_recon', 'para')
        else
            save(sprintf('%s%s_ttv_%.2g_stv_%.2g_frametime_%.2g.mat', path_yt_recon, all_mat(ifile).name(1:end-4), para.weight_tTV, para.weight_sTV, kSpace_info.user_TR * n_tr_per_frame_per_rr / 1000), 'Image_recon', 'para')
        end
    end
end
