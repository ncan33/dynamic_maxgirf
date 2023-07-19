function recon_2D_ungated_SPGR_mSMS_low_dose_KWIC(kSpace_all,kSpace_info,para)

%% recon PD
%kSpace_all = yExportKspace(kSpace_all,0,0.5);
set = kSpace_info.set==0;
set = find(set); set = set(1:600);
nset = max(kSpace_info.set(:))+1;
% RING trajectory correction using PD rays
para.trajectory_correction = RING_SMS(kSpace_all(:,set,:),kSpace_info.angle_mod(set));

if kSpace_info.NumberOfPDlines~=0
    kSpace_info.NumberOfPDlines = kSpace_info.Protocol.lProtonDensMap*nset*kSpace_info.Protocol.sKSpace.lRadialViews;
    para.nor_sl = para.nor_sl;
    Image_PD = recon_PD_multiple_frames(kSpace_all(:,1:kSpace_info.NumberOfPDlines,:),kSpace_info,para);
    save(fullfile(para.dir.save_recon_img_mat_dir,para.dir.save_recon_img_name),'Image_PD','-v7.3')
else
    set = kSpace_info.set==0;
    set = find(set); set = set(361:960);
    % RING trajectory correction using steady_state rays
    para.trajectory_correction = RING_SMS(kSpace_all(:,set,:),kSpace_info.angle_mod(set));
end

%% cut PD rays and non_steaty_state rays
non_steady_state_rays = 0;
kSpace_all(:,1:kSpace_info.NumberOfPDlines+non_steady_state_rays,:) = [];
kSpace_info.angle_mod(:,1:kSpace_info.NumberOfPDlines+non_steady_state_rays,:) = [];
kSpace_info.phase_mod(:,1:kSpace_info.NumberOfPDlines+non_steady_state_rays,:) = [];
kSpace_info.set(:,1:kSpace_info.NumberOfPDlines+non_steady_state_rays,:) = [];
kSpace_info.TimeStamp(:,1:kSpace_info.NumberOfPDlines+non_steady_state_rays,:) = [];
para.kSpace_info = kSpace_info;

[sx,nor_all,no_comp] = size(kSpace_all);

para.Recon.nSMS = max(kSpace_info.phase_mod(:))+1;
kCenter = para.kSpace_center;
%interp_method       = para.Recon.interp_method;
para.setting.ifplot = 0;



%% pre-interp data

for i=1:nset
    set = kSpace_info.set==i-1;
    kSpace_radial = kSpace_all(:,set,:);
    theta = kSpace_info.angle_mod(set);
    phase = kSpace_info.phase_mod(set);
    
    nor_sl = para.nor_sl;
    nor_total = size(kSpace_radial,2);
    nof = floor(nor_total/nor_sl);
    nor_total = nof*nor_sl;
    kSpace_radial(:,nor_total+1:end,:) = [];
    theta(nor_total+1:end) = [];
    phase(nor_total+1:end) = [];
    
    kSpace_radial = reshape(kSpace_radial,[sx,nor_sl,nof,no_comp]);
    theta = reshape(theta,[1,nor_sl,nof]);
    phase = reshape(phase,[1,nor_sl,nof]);
    
    [kx,ky] = get_k_coor(sx,theta,0,kCenter);

    correction = para.trajectory_correction.*permute([cos(theta);sin(theta)],[4,1,2,3]);
    correction = squeeze(sum(correction,2));
    kx = kx - correction(1,:,:);
    ky = ky - correction(2,:,:);
    %phase = get_phase(phase);
    
    [Data{i}.kSpace,Data{i}.G] = GROG.SMS_GROG_KWIC(kSpace_radial,kx,ky,phase,para);
    para.Recon.nor = nor_sl;
    [Data{i},para] = get_Data_SMS(Data{i},para);

    para.Recon.no_comp = no_comp;
end

[para.kSpace_info.TimeStamp] = para.kSpace_info.TimeStamp(1:nset*nor_sl:end);

siz = size(Data{1}.first_est);
siz(4) = nset;
Image = zeros(siz,'single');


%% recon reference image
para.setting.ifGPU = 1;
para.setting.ifplot = 0;
para.Recon.type = 'seperate SMS test';
para.Recon.noi = 70;

for i=1:nset
    scale = max(abs(Data{i}.first_est(:)));
    para.Recon.weight_tTV = scale*para.weight_tTV;
    para.Recon.weight_sTV = scale*para.weight_sTV;
    for j=1:para.Recon.no_comp
        k_temp = Data{i}.kSpace(:,:,:,j,:,:,:);
        kSpace(:,j) = k_temp(Data{i}.mask);
    end
    Data{i}.kSpace = kSpace;
    Image(:,:,:,i,:) = STCR_conjugate_gradient_low_rank_bins(Data{i},para);
    clear kSpace
end

order = vec((1:nset)'+([1,para.Recon.nSMS:-1:2]-1)*nset);
% Image_sys = Image_sys(:,:,:,order);
Image = Image(:,:,:,order);
Image = crop_half_FOV(abs(Image));
%Image = Image(:,:,:,[1,2,3,7,8,9,4,5,6]);
if isfile(fullfile(para.dir.save_recon_img_mat_dir,para.dir.save_recon_img_name))
    save(fullfile(para.dir.save_recon_img_mat_dir,para.dir.save_recon_img_name),'Image','para','-append')
else
    save(fullfile(para.dir.save_recon_img_mat_dir,para.dir.save_recon_img_name),'Image','para','-v7.3')
end
return
keyboard
nof_sys = nof;


%% single set
for i=1:para.Recon.nSMS
    cardiac_signal(:,i) = compare_curve_same_image(crop_half_FOV(Image(:,:,:,i)));
end
cardiac_signal = sum(cardiac_signal,2);
sys = local_max(-cardiac_signal);
figure
plot(cardiac_signal)
hold on
plot(sys,cardiac_signal(sys),'o')

para.kSpace_info.TimeStamp_sys = TimeStamp;

para.Recon.noi = 70;
for i=1:nset
    
    set = kSpace_info.set==i-1;
    kSpace_radial = kSpace_all(:,set,:);

    theta = kSpace_info.angle_mod(set);
    phase = kSpace_info.phase_mod(set);
    nor_sl = para.nor_sl;
    nor_total = size(kSpace_radial,2);
    nof = floor(nor_total/nor_sl);
    nor_total = nof*nor_sl;
    kSpace_radial(:,nor_total+1:end,:) = [];
    theta(nor_total+1:end) = [];
    phase(nor_total+1:end) = [];
    
    kSpace_radial = reshape(kSpace_radial,[sx,nor_sl,nof,no_comp]);
    theta = reshape(theta,[1,nor_sl,nof]);
    phase = reshape(phase,[1,nor_sl,nof]);
    
    [kx,ky] = get_k_coor(sx,theta,0,kCenter);
    correction = para.trajectory_correction.*permute([cos(theta);sin(theta)],[4,1,2,3]);
    correction = squeeze(sum(correction,2));
    kx = kx - correction(1,:,:);
    ky = ky - correction(2,:,:);
    
    Data{i}.kSpace = GROG.SMS_GROG(kSpace_radial,kx,ky,phase,para);
    Data{i}.kSpace = Data{i}.kSpace(:,:,sys(:,i),:,:,:,:);
    Data{i} = get_Data_SMS(Data{i},para);
    scale = max(abs(Data{i}.first_est(:)));
    para.Recon.weight_tTV = scale*para.weight_tTV;
    para.Recon.weight_sTV = scale*para.weight_sTV;
    for j=1:para.Recon.no_comp
        k_temp = Data{i}.kSpace(:,:,:,j,:,:,:);
        kSpace(:,j) = k_temp(Data{i}.mask);
    end
    Data{i}.kSpace = kSpace;
    Image_sys(:,:,:,i,:) = STCR_conjugate_gradient_low_rank_bins(Data{i},para);
    clear kSpace
end
Image_sys = crop_half_FOV(abs(Image_sys));
Image_sys = Image_sys(:,:,:,:);
Image = crop_half_FOV(abs(Image));
Image = Image(:,:,:,:);
order = vec((1:nset)'+([1,para.Recon.nSMS:-1:2]-1)*nset);
Image_sys = Image_sys(:,:,:,order);
Image = Image(:,:,:,order);



%% 

for i=1:nset
    cardiac_signal{i} = compare_curve_same_image(crop_half_FOV(Image(:,:,:,i,3)));
    sys{i} = local_max(-cardiac_signal{i});
    nof_sys = min(length(sys{1}),nof);
    para.kSpace_info.TimeStamp_sys{i} = para.kSpace_info.TimeStamp(sys{i});
end

cardiac_signal_temp = zeros(nof_sys,nset);
sys_temp = zeros(nof_sys,nset);
TimeStamp = zeros(nof_sys,nset);
for i=1:nset
    cardiac_signal_temp(:,i) = cardiac_signal{i}(1:nof_sys);
    sys_temp(:,i) = sys{i}(1:nof_sys);
    TimeStamp(:,i) = para.kSpace_info.TimeStamp_sys{i}(1:nof_sys);
end
para.kSpace_info.TimeStamp_sys = TimeStamp;
sys = sys_temp;

para.Recon.noi = 150;
for i=1:nset
    
    set = kSpace_info.set==i-1;
    kSpace_radial = kSpace_all(:,set,:);

    theta = kSpace_info.angle_mod(set);
    phase = kSpace_info.phase_mod(set);
    nor_sl = para.nor_sl;
    nor_total = size(kSpace_radial,2);
    nof = floor(nor_total/nor_sl);
    nor_total = nof*nor_sl;
    kSpace_radial(:,nor_total+1:end,:) = [];
    theta(nor_total+1:end) = [];
    phase(nor_total+1:end) = [];
    
    kSpace_radial = reshape(kSpace_radial,[sx,nor_sl,nof,no_comp]);
    theta = reshape(theta,[1,nor_sl,nof]);
    phase = reshape(phase,[1,nor_sl,nof]);
    
    [kx,ky] = get_k_coor(sx,theta,0,kCenter);
    correction = para.trajectory_correction.*permute([cos(theta);sin(theta)],[4,1,2,3]);
    correction = squeeze(sum(correction,2));
    kx = kx - correction(1,:,:);
    ky = ky - correction(2,:,:);
    
    Data{i}.kSpace = GROG.SMS_GROG(kSpace_radial,kx,ky,phase,para);
    Data{i}.kSpace = Data{i}.kSpace(:,:,sys(:,i),:,:,:,:);
    Data{i} = get_Data_SMS(Data{i},para);
    scale = max(abs(Data{i}.first_est(:)));
    para.Recon.weight_tTV = scale*para.weight_tTV;
    para.Recon.weight_sTV = scale*para.weight_sTV;
    for j=1:para.Recon.no_comp
        k_temp = Data{i}.kSpace(:,:,:,j,:,:,:);
        kSpace(:,j) = k_temp(Data{i}.mask);
    end
    Data{i}.kSpace = kSpace;
    Image_sys(:,:,:,i,:) = STCR_conjugate_gradient_low_rank_bins(Data{i},para);
    clear kSpace
end
Image_sys = crop_half_FOV(abs(Image_sys));
Image_sys = Image_sys(:,:,:,:);
Image = crop_half_FOV(abs(Image));
Image = Image(:,:,:,:);
order = vec((1:nset)'+([1,para.Recon.nSMS:-1:2]-1)*nset);
Image_sys = Image_sys(:,:,:,order);
Image = Image(:,:,:,order);
%Image = Image(:,:,:,[1,2,3,7,8,9,4,5,6]);
if isfile(fullfile(para.dir.save_recon_img_mat_dir,para.dir.save_recon_img_name))
    save(fullfile(para.dir.save_recon_img_mat_dir,para.dir.save_recon_img_name),'Image_sys','Image','para','-append')
else
    save(fullfile(para.dir.save_recon_img_mat_dir,para.dir.save_recon_img_name),'Image_sys','Image','para','-v7.3')
end
return
