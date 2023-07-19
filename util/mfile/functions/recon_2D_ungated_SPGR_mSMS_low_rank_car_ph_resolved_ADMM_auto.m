function recon_2D_ungated_SPGR_mSMS_low_rank_car_ph_resolved_ADMM_auto(kSpace_all,kSpace_info,para)
tic
%% recon PD
%kSpace_all = yExportKspace(kSpace_all,0,0.5);
save(fullfile(para.dir.save_recon_img_mat_dir,para.dir.save_recon_img_name),'para','-v7.3')

set = kSpace_info.set==0;
set = find(set); set = set(1:600);
nset = max(kSpace_info.set(:))+1;
%% RING trajectory correction using PD rays
para.trajectory_correction = RING_SMS(kSpace_all(:,set,:),kSpace_info.angle_mod(set));

if kSpace_info.NumberOfPDlines~=0
    kSpace_info.NumberOfPDlines = kSpace_info.Protocol.lProtonDensMap*nset*kSpace_info.Protocol.sKSpace.lRadialViews;
%     para.nor_sl = para.nor_sl;
     Image_PD = recon_PD_multiple_frames(kSpace_all(:,1:kSpace_info.NumberOfPDlines,:),kSpace_info,para);
     save(fullfile(para.dir.save_recon_img_mat_dir,para.dir.save_recon_img_name),'Image_PD','-v7.3')
else
    set = kSpace_info.set==0;
    set = find(set); set = set(361:960);
    % RING trajectory correction using steady_state rays
    para.trajectory_correction = RING_SMS(kSpace_all(:,set,:),kSpace_info.angle_mod(set));
end

%% cut PD rays and non_steaty_state rays
non_steady_state_rays = 360;
kSpace_all(:,1:kSpace_info.NumberOfPDlines+non_steady_state_rays,:) = [];
kSpace_info.angle_mod(:,1:kSpace_info.NumberOfPDlines+non_steady_state_rays,:) = [];
kSpace_info.phase_mod(:,1:kSpace_info.NumberOfPDlines+non_steady_state_rays,:) = [];
kSpace_info.set(:,1:kSpace_info.NumberOfPDlines+non_steady_state_rays,:) = [];
kSpace_info.TimeStamp(:,1:kSpace_info.NumberOfPDlines+non_steady_state_rays,:) = [];
para.kSpace_info = kSpace_info;

[sx,nor_all,no_comp] = size(kSpace_all);

nSMS = max(kSpace_info.phase_mod(:))+1;
para.Recon.nSMS = nSMS;
para.setting.ifplot = 0;

%% estimate GROG operator 
for iset=1:nset
    set_idx = kSpace_info.set==iset-1;
    kSpace_temp = kSpace_all(:,set_idx,:);
    theta_temp = kSpace_info.angle_mod(set_idx);
    phase_temp = kSpace_info.phase_mod(set_idx);
    for isms = 1:nSMS
        phase_temp_SMS = phase_temp == nSMS-1;
        kSpace_temp_SMS = kSpace_temp(:,phase_temp_SMS,:);
        kSpace_temp_SMS = permute(kSpace_temp_SMS,[1,2,4,3]);
        theta_temp_SMS = theta_temp(phase_temp_SMS);
        [kx,ky] = get_k_coor(sx,theta_temp_SMS,0,sx/2+1);
        [G{iset,isms}.Gx,G{iset,isms}.Gy] = GROG.get_Gx_Gy(kSpace_temp_SMS,kx,ky);
    end
end
para.G = G;
keyboard
[cardiac_signal,resp_signal,para] = sliding_window_low_res_recon_for_self_gating(kSpace_all,para);
para.resp_signal = resp_signal;

%% pre-interp data
% para.nor_sl = para.nor_sl/2;

sys = local_max(-cardiac_signal);
dia = local_max(cardiac_signal);
dia(dia<sys(2)) = [];
dia(dia>sys(end-1)) = [];

para.Recon.self_gating.sys = sys;
para.Recon.self_gating.dia = dia;

para.Recon.self_gating.cardiac_signal = cardiac_signal;
para.Recon.self_gating.sys = sys;
figure
plot(cardiac_signal)
hold on
plot(sys,cardiac_signal(sys),'o')

nor_sl = 6;
nor_one_frame = para.nor_sl;
para.Recon.nor = nor_one_frame;
ray_idx_start_sys = (sys-1)*nor_sl+1;
% ray_idx_end_sys = (sys-1)*nor_sl+nor_one_frame;

ray_idx_start_dia = (dia-1)*nor_sl+1;

Ns2d = round(mean((ray_idx_start_dia - ray_idx_start_sys(2:end-2))/nor_one_frame));
Nd2s = round(mean((ray_idx_start_sys(3:end-1) - ray_idx_start_dia)/nor_one_frame));
Nphase = Ns2d + Nd2s;

Ncycle = length(ray_idx_start_dia);
% Nphase = round(mean((ray_idx_start_sys(3:end-1) - ray_idx_start_sys(2:end-2))/nor_one_frame));

for i=1:nset
    set_temp = kSpace_info.set==i-1;
    kSpace_temp = kSpace_all(:,set_temp,:);
    theta_temp = kSpace_info.angle_mod(set_temp);
    phase_temp = kSpace_info.phase_mod(set_temp);
    Image_temp = squeeze(para.Image(:,:,:,i,:));
    
    kSpace_reorder = zeros(sx,nor_one_frame,Nphase,Ncycle,no_comp);
    theta_reorder = zeros(1,nor_one_frame,Nphase,Ncycle);
    phase_reorder = repmat([0:nSMS-1],[1,nor_one_frame/nSMS*2,Nphase,Ncycle]);
%     phase_reorder = zeros(1,nor_one_frame,Nphase,Ncycle);
    Image_interp = zeros(sx,sx,Nphase,Ncycle,nSMS);
    resp_signal_reorder = zeros(Nphase,Ncycle);
    
    for icycle=1:Ncycle
        sys_start_temp = ray_idx_start_sys(icycle+1);
        dia_start_temp = ray_idx_start_dia(icycle);
        sys_start_next = ray_idx_start_sys(icycle+2);
        
        % put the systole rays
        kSpace_reorder(:,1:nor_one_frame,1,icycle,:) = kSpace_temp(:,sys_start_temp:sys_start_temp+nor_one_frame-1,:);
        theta_reorder(1,1:nor_one_frame,1,icycle) = theta_temp(sys_start_temp:sys_start_temp+nor_one_frame-1);
        phase_reorder(1,1:nor_one_frame,1,icycle) = phase_temp(sys_start_temp:sys_start_temp+nor_one_frame-1);
        
        % put the systole to diastole rays
        nor_in_between = dia_start_temp - sys_start_temp - nor_one_frame + 1;
        n_frame_in_between = Ns2d-1;
        nor_in_between_one_frame = ceil(nor_in_between/n_frame_in_between/nSMS)*nSMS;
        
        iphase = 2;
        for iframe_in_between = 1:n_frame_in_between
            ray_start_temp = sys_start_temp + nor_one_frame + (iframe_in_between - 1)*nor_in_between_one_frame;
            kSpace_reorder(:,1:nor_in_between_one_frame,iphase,icycle,:) = kSpace_temp(:,ray_start_temp:ray_start_temp+nor_in_between_one_frame-1,:);
            theta_reorder(1,1:nor_in_between_one_frame,iphase,icycle) = theta_temp(ray_start_temp:ray_start_temp+nor_in_between_one_frame-1);
            phase_reorder(1,1:nor_in_between_one_frame,iphase,icycle) = phase_temp(ray_start_temp:ray_start_temp+nor_in_between_one_frame-1);
            iphase = iphase + 1;
        end
        
        
        % interp image
        frame_temp = sys(icycle+1):dia(icycle);
        Image_temp_for_interp = Image_temp(:,:,frame_temp,:);
        Image_temp_for_interp = permute(Image_temp_for_interp,[3,1,2,4]);
        Image_temp_for_interp = Image_temp_for_interp(:,:);
        Image_temp_for_interp = interp1(1:length(frame_temp),Image_temp_for_interp,1:(length(frame_temp)-1)/(Ns2d-1):length(frame_temp));
        Image_temp_for_interp = reshape(Image_temp_for_interp,[Ns2d,sx,sx,nSMS]);
        Image_temp_for_interp = permute(Image_temp_for_interp,[2,3,1,4]);
        Image_interp(:,:,1:Ns2d,icycle,:) = Image_temp_for_interp;
        
        % interp resp_signal
        resp_signal_temp = resp_signal(frame_temp);
        resp_signal_reorder(1:Ns2d,icycle) = interp1(1:length(frame_temp),resp_signal_temp,1:(length(frame_temp)-1)/(Ns2d-1):length(frame_temp));
        
        % put the diastole rays
        kSpace_reorder(:,1:nor_one_frame,iphase,icycle,:) = kSpace_temp(:,dia_start_temp:dia_start_temp+nor_one_frame-1,:);
        theta_reorder(1,1:nor_one_frame,iphase,icycle) = theta_temp(dia_start_temp:dia_start_temp+nor_one_frame-1);
        phase_reorder(1,1:nor_one_frame,iphase,icycle) = phase_temp(dia_start_temp:dia_start_temp+nor_one_frame-1);
        iphase = iphase + 1;
        
        % put the diastole to systole rays
        nor_in_between = sys_start_next - dia_start_temp - nor_one_frame + 1;
        n_frame_in_between = Nd2s-1;
        nor_in_between_one_frame = ceil(nor_in_between/n_frame_in_between/nSMS)*nSMS;
        
        for iframe_in_between = 1:n_frame_in_between
            ray_start_temp = dia_start_temp + nor_one_frame + (iframe_in_between - 1)*nor_in_between_one_frame;
            kSpace_reorder(:,1:nor_in_between_one_frame,iphase,icycle,:) = kSpace_temp(:,ray_start_temp:ray_start_temp+nor_in_between_one_frame-1,:);
            theta_reorder(1,1:nor_in_between_one_frame,iphase,icycle) = theta_temp(ray_start_temp:ray_start_temp+nor_in_between_one_frame-1);
            phase_reorder(1,1:nor_in_between_one_frame,iphase,icycle) = phase_temp(ray_start_temp:ray_start_temp+nor_in_between_one_frame-1);
            iphase = iphase + 1;
        end
        
        % interp image
        frame_temp = dia(icycle)+1:sys(icycle+2)-1;
        Image_temp_for_interp = Image_temp(:,:,frame_temp,:);
        Image_temp_for_interp = permute(Image_temp_for_interp,[3,1,2,4]);
        Image_temp_for_interp = Image_temp_for_interp(:,:);
        Image_temp_for_interp = interp1(1:length(frame_temp),Image_temp_for_interp,1:(length(frame_temp)-1)/(Nd2s-1):length(frame_temp));
        Image_temp_for_interp = reshape(Image_temp_for_interp,[Nd2s,sx,sx,nSMS]);
        Image_temp_for_interp = permute(Image_temp_for_interp,[2,3,1,4]);
        Image_interp(:,:,Ns2d+1:Nphase,icycle,:) = Image_temp_for_interp;
        
        % interp resp_signal
        resp_signal_temp = resp_signal(frame_temp);
        resp_signal_reorder(Ns2d+1:Nphase,icycle) = interp1(1:length(frame_temp),resp_signal_temp,1:(length(frame_temp)-1)/(Ns2d-1):length(frame_temp));
        
    end
    
    % pre-interpolate
    nor_max = size(kSpace_reorder,2);
    kSpace_reorder = reshape(kSpace_reorder,[sx,nor_max,Nphase*Ncycle,no_comp]);
    theta_reorder = reshape(theta_reorder,[1,nor_max,Nphase*Ncycle]);
    phase_reorder(:,nor_max+1:end,:,:) = [];
    phase_reorder = reshape(phase_reorder,[1,nor_max,Nphase*Ncycle]);
    
    [kx,ky] = get_k_coor(sx,theta_reorder,0,sx/2+1);
    
    % RING trajectory correction
    correction = para.trajectory_correction.*permute([cos(theta_reorder);sin(theta_reorder)],[4,1,2,3]);
    correction = squeeze(sum(correction,2));
    kx = kx - correction(1,:,:);
    ky = ky - correction(2,:,:);
    
    for isms=1:nSMS
        G{isms} = para.G{i,isms};
    end
    Data{i}.kSpace = GROG.SMS_GROG(kSpace_reorder,kx,ky,phase_reorder,para,G);
    Data{i}.first_guess = reshape(Image_interp,[sx,sx,Nphase*Ncycle,1,nSMS]);
    [Data{i},para] = get_Data_SMS(Data{i},para);
    
end
    
para.Recon.self_gating.resp_signal = resp_signal_reorder(:);

para.Recon.no_comp = no_comp;
para.Recon.bins = repmat(diag(true(1,Nphase)),[1,Ncycle]);
% im_bin = reshape(im_bin,[sx,sx,nof,nset,nSMS]);

para.setting.ifGPU = 0;
para.setting.ifplot = 0;
para.Recon.noi = 50;
para = rmfield(para,'Image');
save(fullfile(para.dir.save_recon_img_mat_dir,para.dir.save_recon_img_name),'para','-append')

%% 2nd reconstruction stage

% Data_all = rmfield(Data{1},{'G','first_est'});
Data_all = Data{1};
scale = max(abs(Data{1}.first_est(:)));
for i=2:nset
    Data_all.kSpace = cat(6,Data_all.kSpace,Data{i}.kSpace);
    Data_all.mask = cat(6,Data_all.mask,Data{i}.mask);
    scale = max([scale;abs(Data{i}.first_est(:))]);
    Data_all.sens_map = cat(6,Data_all.sens_map,Data{i}.sens_map);
    Data_all.first_guess = cat(6,Data_all.first_guess,Data{i}.first_guess);
end
clear Data
for i=1:para.Recon.no_comp
    k_temp = Data_all.kSpace(:,:,:,i,:,:,:);
    kSpace(:,i) = k_temp(Data_all.mask);
end
Data_all.kSpace = kSpace;
clear kSpace k_temp

% Data_all.first_guess = permute(im_bin,[1,2,3,6,5,4]);
para.Recon.weight_tTV = scale*0.04;
para.Recon.weight_sTV = scale*0.00;
para.Recon.type = 'seperate SMS test';
para.Recon.noi = 50;

% Data_all.myo_mask = mask;
clearvars -except Data_all para

%[Data_all,para] = get_iso_image(Data_all,para);
%Data_all = Patch_tracking(Data_all,para);

para.setting.ifGPU = 1;
% Image = STCR_conjugate_gradient(Data_all,para);
% Data_all.first_guess = permute(Image,[1,2,3,6,4,5]);

Data_all = rmfield(Data_all,'first_est');
[sx,sy,nof,~,nSMS,nset] = size(Data_all.first_guess);
nslice = nSMS*nset;
order = 1:nSMS:nslice;
for i=nSMS:-1:2
    order = [order,i:nSMS:nslice];
end
[~,order_back] = sort(order);

Data_all.first_guess = reshape(Data_all.first_guess,[sx,sy,nof,nslice]);
Data_all.first_guess = Data_all.first_guess(:,:,:,order);
Data_all.first_guess = single(Data_all.first_guess);

dia_loc = round(size(para.Recon.bins,1)/2+1);
Image_dia = abs(crop_half_FOV(Data_all.first_guess(:,:,para.Recon.bins(dia_loc,:),:)));
Image_dia = permute(Image_dia,[1,2,4,3]);

%
para.mask = imresize(para.mask,[sx/2,sx/2]);
para.mask = para.mask(:,:,order);

[~,shifts_dia,~] = rigid_reg_GS_3D(Image_dia,para.mask);



% resp_signal = zeros(size(Data_all.first_guess,3),1);
% for i=1:size(Data_all.first_guess,4)
%     resp_signal = resp_signal + compare_curve_same_image(Data_all.first_guess(:,:,:,i));
% end
% resp_signal = smooth(resp_signal,15);
% 
% resp_signal_correction = polyfit(1:length(resp_signal),resp_signal',1);
% resp_signal_correction = polyval(resp_signal_correction,1:length(resp_signal));
% resp_signal = resp_signal./resp_signal_correction';

figure
plot(para.resp_signal)

shifts_interp = interp_rigid_shifts(shifts_dia,para.Recon.self_gating.resp_signal,para.Recon.bins(dia_loc,:));
shifts_interp = shifts_interp - mean(shifts_interp,2);


Data_all.shifts = shifts_interp;
para.shifts = Data_all.shifts;
% para.myo_mask = Data_all.myo_mask;

% mask = cat(1,false(sx/2,sy/2,nslice),mask);
% mask = cat(2,false(sx,sy/2,nslice),mask);
% mask = circshift(mask,[-sx/4,-sx/4]);
% 
% Data_all.myo_mask = mask;
% Data_all.llr = Patch_init_with_guide(Data_all,para,mask,shifts_interp,[5,5,5],3);
% Data_all.llr_out = Patch_init_with_guide(Data_all,para,~mask,zeros(size(shifts_interp)),[5,5,5],5);

Data_all.iso.Nphase = size(para.Recon.bins,1);
Data_all.iso.Ncycle = nof/Data_all.iso.Nphase;
Data_all.iso.order = order;
Data_all.iso.order_back = order_back;
Data_all.iso.nSMS = nSMS;
Data_all.iso.nset = nset;


clearvars -except Data_all para

% 
% Data_all.first_guess = interp_2_slices_in_between(permute(Data_all.first_guess,[1,2,4,3]));
% Data_all.first_guess = permute(Data_all.first_guess,[1,2,4,3]);

% Data_all.myo_mask = interp_2_slices_in_between(Data_all.myo_mask);
% Data_all.myo_mask = logical(Data_all.myo_mask);
% Data_all.llr = Patch_tracking_4D_with_guide(Data_all,para,Data_all.myo_mask,Data_all.shifts);

%% 4D
% patch_size = [5,5,2,9];
% search_size = [9,9,4,9];
% patch_shift = [5,5,2];
% Data_all.llr = Patch_tracking_4D_with_guide_random_search(Data_all.first_guess,Data_all.shifts,patch_size,search_size,patch_shift);
% 
% para.setting.ifGPU = 1;
% para.Recon.noi = 50;
% 
% im_bin = LLR_patch_tracking_4D_ADMM(Data_all,para);
%%

%% 3D
patch_size = [5,5,2];
search_size = [9,9,4];
patch_shift = [5,5,2];
Nphase = size(para.Recon.bins,1);
Data_all.llr = Patch_tracking_3D_with_guide_random_search(Data_all.first_guess,Nphase,Data_all.shifts,patch_size,search_size,patch_shift);

para.setting.ifGPU = 1;
para.Recon.noi = 50;
toc;
keyboard
im_bin = LLR_patch_tracking_4D_ADMM(Data_all,para);
%%
% im_bin = STCR_LLR_patch_tracking_POCS(Data_all,para);
% 
% 
% im_bin = LLR_patch_tracking_3D_ADMM(Data_all,para);
% 
% para.Recon.noi = 50;
% im_bin = STCR_conjugate_gradient_MSMS_iso_patch_llrt(Data_all,para);


% 
% im_bin = STCR_conjugate_gradient_MSMS_iso_patch_llr(Data_all,para);
im_bin = permute(im_bin,[1,2,4,3]);
im_bin = abs(crop_half_FOV(im_bin));
save(fullfile(para.dir.save_recon_img_mat_dir,para.dir.save_recon_img_name),'im_bin','para','-append')

Image_sys = im_bin(:,:,:,para.Recon.bins(1,:));
% Image_sys = permute(Image_sys,[1,2,4,3]);
save(fullfile(para.dir.save_recon_img_mat_dir,para.dir.save_recon_img_name),'Image_sys','-append')

dia_loc = round((size(para.Recon.bins,1)+1)/2);
Image_dia = im_bin(:,:,:,para.Recon.bins(dia_loc,:));
% Image_dia = permute(Image_dia,[1,2,4,3]);
save(fullfile(para.dir.save_recon_img_mat_dir,para.dir.save_recon_img_name),'Image_dia','-append')

return






      