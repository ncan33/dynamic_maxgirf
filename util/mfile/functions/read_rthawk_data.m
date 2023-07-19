function read_rthawk_data(area, vol_number, file_index)
if nargin == 2
    file_index = [];
end

if ismac
    path_server = dir(sprintf('/Users/ytian/MRIdata/%s/*%d*', area, vol_number));
    path_server = fullfile(path_server.folder, path_server.name);
    path_yt = path_server;
else
    if contains(area, 'phantom')
        path_server = dir(sprintf('/server/sdata/ytian/mri_data/disc/phantom/*%d*', vol_number));
    else
        path_server = dir(sprintf('/server/sdata/mri_data/disc/%s/vol%04g_*', area, vol_number));
    end
    path_server = fullfile(path_server.folder, path_server.name);
    
    path_yt = strfind(path_server, 'mri_data');
    path_yt = path_server(path_yt:end);
    path_yt = ['/server/sdata/ytian/', path_yt];
end


if ~isfolder(path_yt)
    mkdir(path_yt)
end

path_yt_raw = [path_yt, '/raw_hawk/'];

if ~isfolder(path_yt_raw)
    mkdir(path_yt_raw)
end

path_server_raw = [path_server, '/raw_hawk/'];

path_search_all = {};
path_search_all{1} = path_server_raw;
sub_folder_temp = search_sub_folder(path_search_all{1});
sub_folder_to_be_searched = sub_folder_temp;
sub_folder_all = [path_server_raw; sub_folder_temp];
while ~isempty(sub_folder_to_be_searched)
    sub_folder_temp = search_sub_folder(sub_folder_to_be_searched{1});
    sub_folder_to_be_searched(1) = [];
    sub_folder_to_be_searched = [sub_folder_to_be_searched; sub_folder_temp];
    sub_folder_all = [sub_folder_all; sub_folder_temp];
end

all_dat = [];
for K = 1 : length(sub_folder_all)
    thisdir = sub_folder_all{K};
    dat_temp = dir(fullfile(thisdir, 'usc_disc_yt*raw.dat'));
    all_dat = [all_dat; dat_temp];
end

% all_dat = dir([path_server_raw, '*raw.dat']);

n_files = length(all_dat);
if isempty(file_index)
    file_index = 1:n_files;
end

for i = file_index %1:n_files
    
    try
        t1 = tic;
        file_name = fullfile(all_dat(i).folder, all_dat(i).name);
        fprintf('--------------------------------- \n')
        fprintf('reading data:  %d/%d \n', i, n_files)
        fprintf('file name:     %s \n',file_name)
        
        [kspace, header] = loadRthData3(file_name);
        
        kspace = single(kspace);
        
        %% trajectory
        if isfield(header(1), 'kspace')
            ksDims = header(1).kspace.dimensions;
        else
            ksDims = 2;
        end
        
        nSamples = header(1).extent(1);
        nCoils = header(1).extent(2);
        nArmsTotal  = size(kspace, 3); %length(header);
        nSamples = size(kspace, 2) / nCoils;
        
        %% kspace data
        kspace = squeeze(kspace(1, :, :) + 1i * kspace(2, :, :));
        kspace = reshape(kspace, [nSamples, nCoils, nArmsTotal]);
        
        %% drop some data at end
        %     nFrames = floor(nArmsTotal / GAsteps);
        %     data(:, :, GAsteps * nFrames + 1 : end) = [];
        
        %% load uncorrected trajectory
        file_name = fullfile(all_dat(i).folder, [all_dat(i).name(1:end-7), 'traj.dat']);
        [~, header_traj, trajectory] = loadRthData3(file_name);
        trajectory = single(trajectory);
        header_traj = header_traj(1);
        
        GAsteps = size(trajectory, 2) / nSamples / (ksDims + 1);
        ntraj   = size(trajectory, 3); % number of trajectory
        
        %% girf correction
        q0 = header_traj.user_QuaternionW;
        q1 = header_traj.user_QuaternionX;
        q2 = header_traj.user_QuaternionY;
        q3 = header_traj.user_QuaternionZ;
        
        rot = [2 * (q0^2  + q1^2 ) - 1,     2 * (q1*q2 - q0*q3),        2 * (q1*q3 + q0*q2);
            2 * (q1*q2 + q0*q3),         2 * (q0^2  + q2^2 ) - 1,    2 * (q2*q3 - q0*q1);
            2 * (q1*q3 - q0*q2),         2 * (q2*q3 + q0*q1),        2 * (q0^2  + q3^2) - 1];
        
        try
            nSamples_traj = nSamples;
            trajectory = reshape(trajectory, [ksDims + 1, nSamples_traj, GAsteps, ntraj]);
        catch
            GAsteps = header_traj.user_interleaves;
            nSamples_traj = size(trajectory, 2) / GAsteps / (ksDims + 1);
            trajectory = reshape(trajectory, [ksDims + 1, nSamples_traj, GAsteps, ntraj]);
        end
        %     trajectory = permute(trajectory, [1, 2, 4, 3]);
        %     trajectory = reshape(trajectory, [ksDims + 1, nSamples * 2, 2, GAsteps]);
        %     trajectory = permute(trajectory, [1, 2, 4, 3]);
        
        if ntraj == 1
            kx0 = permute(trajectory(1, :, :, 1), [2, 3, 1]);
            ky0 = permute(trajectory(2, :, :, 1), [2, 3, 1]);
            kx  = permute(trajectory(1, :, :, 1), [2, 3, 1]); % corrected traj
            ky  = permute(trajectory(2, :, :, 1), [2, 3, 1]); % corrected traj
            w = permute(trajectory(3, :, :, 1), [2, 3, 1]);
        elseif ntraj == 2
            kx0 = permute(trajectory(1, :, :, 1), [2, 3, 1]);
            ky0 = permute(trajectory(2, :, :, 1), [2, 3, 1]);
            kx  = permute(trajectory(1, :, :, 1), [2, 3, 1]); % corrected traj
            ky  = permute(trajectory(2, :, :, 1), [2, 3, 1]); % corrected traj
            w = permute(trajectory(3, :, :, 1), [2, 3, 1]);
            ntraj = 1;
        elseif ntraj == 3
            kx0 = permute(trajectory(1, :, :, :), [2, 3, 1, 4]);
            ky0 = permute(trajectory(2, :, :, :), [2, 3, 1, 4]);
            kx  = squeeze(kx0);
            ky  = squeeze(ky0);
            kx0 = reshape(kx0, [nSamples, GAsteps * ntraj]);
            ky0 = reshape(ky0, [nSamples, GAsteps * ntraj]);
            w   = permute(trajectory(3, :, :, :), [2, 3, 4, 1]);
        end
        
        
        Gx0 = [kx0(1, :); diff(kx0)];
        Gy0 = [ky0(1, :); diff(ky0)];
        
        %     load('./WaveForm2/waveForm.mat', 'gx', 'gy')
        
        %     Gx0 = gx;
        %     Gy0 = gy;
        
        WaveForm = cat(3, Gx0, Gy0);
        WaveForm(:, :, 3) = 0;
        
        dt = 1 / header_traj.user_samplingRate / 1000;
        
        R.T = 0.55;
        R.R = rot;
        
        %     [kPred, GPred] = apply_GIRF(WaveForm, dt, R);
        kPred = apply_GIRF_20210827(WaveForm, dt, R, -6, 0);
        
        kx_GIRF = kPred(:, :, 1);
        ky_GIRF = kPred(:, :, 2);
        
        if ntraj < 3
            k_scale = max(vec(sqrt(kx_GIRF.^2 + ky_GIRF.^2)));
            kx_GIRF = kx_GIRF ./ k_scale / 2; % * matrix_size(1);
            ky_GIRF = ky_GIRF ./ k_scale / 2; % * matrix_size(2);
        elseif ntraj == 3
            kx_GIRF = reshape(kx_GIRF, [nSamples_traj, GAsteps, ntraj]);
            ky_GIRF = reshape(ky_GIRF, [nSamples_traj, GAsteps, ntraj]);
            k_scale = max(max(sqrt(kx_GIRF.^2 + ky_GIRF.^2)));
            kx_GIRF = kx_GIRF ./ k_scale / 2; % * matrix_size(1);
            ky_GIRF = ky_GIRF ./ k_scale / 2; % * matrix_size(2);
        end
        
        kx_GIRF = single(kx_GIRF);
        ky_GIRF = single(ky_GIRF);
        
        %% demodultion delay phase
    dx = header_traj.user_TranslationX;
    dy = header_traj.user_TranslationY;
    dz = header_traj.user_TranslationZ;
    
    % trjectory to physical unit
    kx0_demod = kx0 / header_traj.user_ResolutionX;
    ky0_demod = ky0 / header_traj.user_ResolutionY;
    
    kx_girf_demod = kx_GIRF / header_traj.user_ResolutionX;
    ky_girf_demod = ky_GIRF / header_traj.user_ResolutionY;
    
    % rotate trajectory to physical coordinate
    k0 = cat(1, kx0_demod(:)', ky0_demod(:)', zeros(size(kx0_demod(:)))');
    k0 = rot * k0;
    
    kx0_demod = reshape(k0(1, :), [nSamples, GAsteps * ntraj]);
    ky0_demod = reshape(k0(2, :), [nSamples, GAsteps * ntraj]);
    kz0_demod = reshape(k0(3, :), [nSamples, GAsteps * ntraj]);
    
    k_girf = cat(1, kx_girf_demod(:)', ky_girf_demod(:)', zeros(size(kx_girf_demod(:)))');
    k_girf = rot * k_girf;
    
    kx_girf_demod = reshape(k_girf(1, :), [nSamples, GAsteps * ntraj]);
    ky_girf_demod = reshape(k_girf(2, :), [nSamples, GAsteps * ntraj]);
    kz_girf_demod = reshape(k_girf(3, :), [nSamples, GAsteps * ntraj]);
    
    
    phase_x_0 = 2 * pi * dx * kx0_demod;
    phase_y_0 = 2 * pi * dy * ky0_demod;
    phase_z_0 = 2 * pi * dz * kz0_demod;
    
    phase_x_girf = 2 * pi * dx * kx_girf_demod;
    phase_y_girf = 2 * pi * dy * ky_girf_demod;
    phase_z_girf = 2 * pi * dz * kz_girf_demod;
    
    demod_phase_x = circshift(phase_x_0, [-2, 0]) - phase_x_girf;
    demod_phase_y = circshift(phase_y_0, [-2, 0]) - phase_y_girf;
    demod_phase_z = circshift(phase_z_0, [-2, 0]) - phase_z_girf;
    
    demod_phase = demod_phase_x + demod_phase_y + demod_phase_z;
%     delay_corr = reshape(demod_phase(:, viewOrder), [nsample, Narms_per_frame, Nframes]);
    delay_corr = exp(-1i * demod_phase);
        
        %% save
        viewOrder = cat(2, header(:).user_view) + 1;
        
        kspace = permute(kspace, [1, 3, 2]);
        kspace = reshape(kspace, [nSamples, nArmsTotal, nCoils]);
        
        kspace_info = header_traj;
        kspace_info.kx = kx;
        kspace_info.ky = ky;
        kspace_info.kx_GIRF = kx_GIRF;
        kspace_info.ky_GIRF = ky_GIRF;
        kspace_info.DCF = w;
        kspace_info.demodultion_delay = delay_corr;
        
        if isfield(header_traj, 'user_nSliceScan') || isfield(header_traj, 'user_Slices')
            if isfield(header_traj, 'user_nSliceScan')
                nslice = header_traj.user_nSliceScan;
                idx_slice = cat(2, header(:).user_RFIndex);
                nslice = max(idx_slice) + 1;
            else
                nslice = header_traj.user_Slices;
                idx_slice = cat(2, header(:).user_Slices);
                nslice = max(idx_slice) + 1;
                if contains(file_name, 'sms')
                    nslice = header_traj.user_Slices / 3;
                    idx_slice = repmat(1:nslice, [nArmsTotal/nslice, 1]);
                    idx_slice = idx_slice(:);
                end
                
            end
            
%             idx_slice = idx_slice - min(idx_slice) + 1;
            idx_slice = idx_slice + 1;
            kspace_all = kspace;
            for islice = 1:nslice
                idx = idx_slice == islice;
                kspace_info.viewOrder = viewOrder(idx);
                kspace_info.timeStamp = cat(1, header(idx).timeStamp);
                
                if isfield(header, 'user_RFIndex')
                    kspace_info.RFIndex = cat(2, header(idx).user_RFIndex);
                end
                if isfield(header, 'user_TimeSinceTrig')
                    kspace_info.TimeSinceTrig = cat(2, header(idx).user_TimeSinceTrig);
                end
                if isfield(header, 'user_viewSegmentIndex')
                    kspace_info.ViewSegmentIndex = cat(2, header(idx).user_viewSegmentIndex);
                end
                if isfield(header, 'user_ResolutionIndex')
                    kspace_info.ResolutionIndex = cat(2, header(idx).user_ResolutionIndex);
                end
                if isfield(header, 'user_QuaternionX')
                    kspace_info.QuaternionX = cat(2, header(idx).user_QuaternionX);
                    kspace_info.QuaternionY = cat(2, header(idx).user_QuaternionY);
                    kspace_info.QuaternionZ = cat(2, header(idx).user_QuaternionZ);
                    kspace_info.QuaternionW = cat(2, header(idx).user_QuaternionW);
                end
                
                raw_dir = file_name;
                kspace = kspace_all(:, idx, :);
                if ~isempty(kspace)
                    save(sprintf('%s%s_slice_%02g.mat', path_yt_raw, all_dat(i).name(1:end-8), islice), 'kspace', 'kspace_info', 'raw_dir',  '-v7.3')
                end
            end
        else
            
            kspace_info.viewOrder = viewOrder;
            kspace_info.timeStamp = cat(1, header(:).timeStamp);
            
            if isfield(header, 'user_RFIndex')
                kspace_info.RFIndex = cat(2, header(:).user_RFIndex);
            end
            if isfield(header, 'user_TimeSinceTrig')
                kspace_info.TimeSinceTrig = cat(2, header(:).user_TimeSinceTrig);
            end
            if isfield(header, 'user_viewSegmentIndex')
                kspace_info.ViewSegmentIndex = cat(2, header(:).user_viewSegmentIndex);
            end
            if isfield(header, 'user_ResolutionIndex')
                kspace_info.ResolutionIndex = cat(2, header(:).user_ResolutionIndex);
            end
            if isfield(header, 'user_QuaternionX')
                kspace_info.QuaternionX = cat(2, header(:).user_QuaternionX);
                kspace_info.QuaternionY = cat(2, header(:).user_QuaternionY);
                kspace_info.QuaternionZ = cat(2, header(:).user_QuaternionZ);
                kspace_info.QuaternionW = cat(2, header(:).user_QuaternionW);
            end
            
            raw_dir = file_name;
            
            save(sprintf('%s%s.mat', path_yt_raw, all_dat(i).name(1:end-8)), 'kspace', 'kspace_info', 'raw_dir', '-v7.3')
        end
        fprintf('done. '), toc(t1), fprintf('\n')
    catch
        fprintf('fail. '), toc(t1), fprintf('\n')
    end
end