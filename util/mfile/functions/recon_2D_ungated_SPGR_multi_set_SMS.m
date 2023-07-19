function [Image_sys,Image_dia,Image,para] = recon_2D_ungated_SPGR_multi_set_SMS(kSpace_all,kSpace_info,para)

kSpace_all(:,1:kSpace_info.NumberOfPDlines,:) = [];
kSpace_info.angle_mod(:,1:kSpace_info.NumberOfPDlines,:) = [];
kSpace_info.phase_mod(:,1:kSpace_info.NumberOfPDlines,:) = [];
kSpace_info.set(:,1:kSpace_info.NumberOfPDlines,:) = [];
kSpace_info.TimeStamp(:,1:kSpace_info.NumberOfPDlines,:) = [];
para.kSpace_info = kSpace_info;

[sx,nor_all,no_comp] = size(kSpace_all);

para.Recon.nSMS = max(kSpace_info.phase_mod(:))+1;
kCenter = para.kSpace_center;
%interp_method       = para.Recon.interp_method;
para.setting.ifplot = 0;

nset = max(kSpace_info.set(:))+1;

%% recon reference 
keyboard
for i=1:nset
    set = kSpace_info.set==i-1;
    kSpace_radial = kSpace_all(:,set,:);
    theta = kSpace_info.angle_mod(set);
    phase = kSpace_info.phase_mod(set);
    
    nor_sl = 9;
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

%     correction = para.trajectory_correction.*permute([cos(theta);sin(theta)],[4,1,2,3]);
%     correction = squeeze(sum(correction,2));
%     kx = kx + correction(1,:,:);
%     ky = ky + correction(1,:,:);
    %phase = get_phase(phase);
    
    [Data{i}.kSpace,Data{i}.G] = GROG.SMS_GROG(kSpace_radial,kx,ky,phase,para);
    para.Recon.nor = nor_sl;
    [Data{i},para] = get_Data_SMS(Data{i},para);

    para.Recon.no_comp = no_comp;
end

para.kSpace_info.TimeStamp = para.kSpace_info.TimeStamp(1:nset*nor_sl:end);

Data_all = Data{1};
for i=2:nset
    Data_all.kSpace = cat(6,Data_all.kSpace,Data{i}.kSpace);
    Data_all.mask = cat(6,Data_all.mask,Data{i}.mask);
    Data_all.first_est = cat(6,Data_all.first_est,Data{i}.first_est);
    Data_all.sens_map = cat(6,Data_all.sens_map,Data{i}.sens_map);
end
clear Data

scale = max(abs(Data_all.first_est(:)));
para.Recon.weight_tTV = scale*0.04;
para.Recon.weight_sTV = scale*0.001;
para.Recon.type = 'seperate SMS new';
para.Recon.noi = 50;

% recon by blocks
siz = size(squeeze(Data_all.first_est));
Image = zeros(siz,'single');

Block_frame = 70;
Nblock = ceil(siz(3)/Block_frame);

if Nblock>1
    Data_temp.sens_map = Data_all.sens_map;
    Data_temp.filter = Data_all.filter;
    Data_temp.SMS = Data_all.SMS;
    for i=1:Nblock
        range = (1:Block_frame) + (i-1)*Block_frame;
        if i==1
            range = [range,max(range)+1:max(range+5)];
        elseif i == Nblock
            range(range>siz(3)) = [];
            range = [min(range)-5:min(range)-1,range];
        else
            range = [min(range)-5:min(range)-1,range,max(range)+1:max(range+5)];
        end
        Data_temp.kSpace = Data_all.kSpace(:,:,range,:,:,:,:);
        Data_temp.mask = Data_all.mask(:,:,range,:,:,:,:);
        Data_temp.first_est = Data_all.first_est(:,:,range,:,:,:,:);
        Image_temp = squeeze(STCR_conjugate_gradient_MSMS(Data_temp,para));
        if i==1
            range(end-4:end) = [];
            flag = 0;
        elseif i == Nblock
            range(1:5) = [];
            flag = 1;
        else
            range(1:5) = [];
            flag = 1;
            range(end-4:end) = [];
        end
        Image(:,:,range,:,:) = Image_temp(:,:,range - (i-1)*Block_frame + flag*5,:,:);
    end
else
    Image = squeeze(STCR_conjugate_gradient_MSMS(Data_all,para));
end
%save([para.dir.save_recon_img_mat_dir,para.dir.save_recon_img_name],'Image','para','-v7.3');
%keyboard
%Image = squeeze(STCR_conjugate_gradient_MSMS(Data_all,para));

%% self gating
%cardiac_signal = compare_curve_same_image(crop_half_FOV(Image(:,:,:,2,2)));

nslice = siz(4)*siz(5);
cardiac_signal = zeros(siz(3),1);
% for i=1:nslice
%     cardiac_signal = cardiac_signal + squeeze(sum(sum(sum(abs(Image(:,:,:,i).*find_LV_RV_2D_ungated_continues(Image(:,:,:,i))))),4));
% end
cardiac_signal = squeeze(sum(sum(sum(sum(abs(crop_half_FOV(Image)))),4),5));
sys = local_max(-cardiac_signal);
dia = local_max(cardiac_signal);
cardiac_signal_smooth = smooth(cardiac_signal,10);
idx_drop = cardiac_signal(sys) > cardiac_signal_smooth(sys);
% sys(idx_drop) = [];
% sys(end) = [];
%idx_drop = cardiac_signal(dia) < cardiac_signal_smooth(dia);
%dia(idx_drop) = [];
dia(1) = [];
figure,plot(cardiac_signal)
hold on
plot(sys,cardiac_signal(sys),'o')
plot(dia,cardiac_signal(dia),'o')
para.cardiac_signal = cardiac_signal;
para.self_gating_systolic = sys;
para.self_gating_diastolic = dia;

RayPosition = sys*nor_sl-round(nor_sl/2);
RayNumber = diff(RayPosition);
RayNumber = [RayPosition(1);RayNumber;nor_total-RayPosition(end)];
idx_sys = false(nof,1);
idx_sys(sys) = true;

%% respiration gating
%{
respiration_signal = compare_curve_same_image(crop_half_FOV(Image(:,:,:,2,3)));
low = local_max(-respiration_signal);
respiration_signal_smooth = smooth(respiration_signal,100);



d_signal = respiration_signal-respiration_signal_smooth;
peaks = local_max(d_signal);
peaks_interp = interp1(peaks,d_signal(peaks),1:nof);
idx = (d_signal)./peaks_interp';
idx = idx(51:150);% blood pool first pass
idx = sort(idx);
idx = idx(100*0.6);
idx = (d_signal)./peaks_interp'>idx;
idx_dia = false(nof,1);
dia = round(para.self_gating_systolic(1:end-1) + (para.self_gating_systolic(2:end) - para.self_gating_systolic(1:end-1))/2);
idx_dia(dia) = true;
cine_cardiac_cycle_pick = idx_dia & idx;
cine_cardiac_cycle_pick = find(cine_cardiac_cycle_pick);
cine_cardiac_cycle_pick(cine_cardiac_cycle_pick<50) = [];
cine_cardiac_cycle_pick(cine_cardiac_cycle_pick>150) = [];
for i=1:length(cine_cardiac_cycle_pick)
    [~,idx] = sort(abs(cine_cardiac_cycle_pick(i) - para.self_gating_systolic));
    cardiac_cycles_begin(i) = idx(1);
    cardiac_cycles_end(i) = idx(2);
end
cine_cardiac_cycle_pick = min(cardiac_cycles_begin,cardiac_cycles_end);
cine_cardiac_cycle_pick = para.self_gating_systolic(cine_cardiac_cycle_pick);

% idx = respiration_signal - respiration_signal_smooth;
% idx_sort = sort(idx);
% idx = idx>idx_sort(round(nof*0.75));
% 
% para.respiratiopn_signal = respiration_signal;
%% cine

% % pick cardiac cycles that have same respiration state
% cine_cardiac_cycle_pick = idx_sys & idx;
% cine_cardiac_cycle_pick = find(cine_cardiac_cycle_pick);
% cine_cardiac_cycle_pick(cine_cardiac_cycle_pick<50) = [];
% cine_cardiac_cycle_pick(cine_cardiac_cycle_pick>150) = [];

number_cine_cycles = length(cine_cardiac_cycle_pick);
kSpace_cine = zeros(sx,1,number_cine_cycles,no_comp,nset);
theta_cine = zeros(1,1,number_cine_cycles,nset);
phase_cine = zeros(1,1,number_cine_cycles,nset);
for i=1:number_cine_cycles
    cycle_begin_ray = cine_cardiac_cycle_pick(i)*nor_sl - round(nor_sl/2);
    cycle_end_ray = sys(find(sys==cine_cardiac_cycle_pick(i))+1)*nor_sl - round(nor_sl/2) - 1;
    selected_rays = cycle_begin_ray:cycle_end_ray;
    nor_temp = floor(length(selected_rays)/nSMS)*nSMS;
    selected_rays(nor_temp+1:end) = [];
    for j=1:nset
        set = kSpace_info.set==j-1;
        kSpace_radial = kSpace_all(:,set,:);
        theta = kSpace_info.angle_mod(set);
        phase = kSpace_info.phase_mod(set);
        kSpace_temp = kSpace_radial(:,selected_rays,:);
        theta_temp = theta(selected_rays);
        phase_temp = phase(selected_rays);
        
        kSpace_cine(:,1:nor_temp,i,:,j) = kSpace_temp;
        theta_cine(1,1:nor_temp,i,j) = theta_temp;
        phase_cine(1,1:nor_temp,i,j) = phase_temp;
    end
end

nor_cine_each_cycle = squeeze(sum(kSpace_cine(sx/2,:,:,1,1)~=0));
nor_cine = floor(median(nor_cine_each_cycle)/nSMS)*nSMS;
cine_window_length = 3;
cine_sliding_length = 3;

nof_cine = (nor_cine-cine_window_length+cine_sliding_length)/cine_sliding_length;

% resort rays from selected cardiac cycles
kSpace_cine_resort = zeros(sx,1,nof_cine,no_comp,nset);
theta_cine_resort = zeros(1,1,nof_cine,nset);
phase_cine_resort = zeros(1,1,nof_cine,nset);
for i=1:nof_cine
    if i==nof_cine
        selected_rays_start = cine_window_length + (i-1)*cine_sliding_length;
        for j=1:nset
            kSpace_temp = kSpace_cine(:,selected_rays_start:end,:,:,j);
            theta_temp = theta_cine(1,selected_rays_start:end,:,j);
            phase_temp = phase_cine(1,selected_rays_start:end,:,j);
            nor_temp = length(phase_temp(:));
            kSpace_temp = reshape(kSpace_temp,[sx,nor_temp,1,no_comp]);
            theta_temp = theta_temp(:).';
            phase_temp = phase_temp(:).';
            kSpace_cine_resort(:,1:nor_temp,i,:,j) = kSpace_temp;
            theta_cine_resort(1,1:nor_temp,i,j) = theta_temp;
            phase_cine_resort(1,1:nor_temp,i,j) = phase_temp;
        end
    else
        selected_rays = (1:cine_window_length) + (i-1)*cine_sliding_length;
        for j=1:nset
            kSpace_temp = kSpace_cine(:,selected_rays,:,:,j);
            theta_temp = theta_cine(1,selected_rays,:,j);
            phase_temp = phase_cine(1,selected_rays,:,j);
            nor_temp = length(phase_temp(:));
            kSpace_temp = reshape(kSpace_temp,[sx,nor_temp,1,no_comp]);
            theta_temp = theta_temp(:).';
            phase_temp = phase_temp(:).';
            kSpace_cine_resort(:,1:nor_temp,i,:,j) = kSpace_temp;
            theta_cine_resort(1,1:nor_temp,i,j) = theta_temp;
            phase_cine_resort(1,1:nor_temp,i,j) = phase_temp;
        end
    end
end
% gridding
for i=1:nset
    kSpace_temp = kSpace_cine_resort(:,:,:,:,i);
    theta_temp = theta_cine_resort(:,:,:,i);
    phase_temp = phase_cine_resort(:,:,:,i);
    idx_no_data = kSpace_temp(sx/2,:,:,1)==0;
    [kx,ky] = get_k_coor(sx,theta_temp,0,kCenter);
    for j=1:nSMS
        SMS_idx = phase_temp==j-1;
        SMS_idx = SMS_idx & ~idx_no_data;
        for iframe = 1:nof_cine
            SMS_idx_temp = SMS_idx(:,:,iframe);
            kSpace_temp_frame = kSpace_temp(:,SMS_idx_temp,iframe,:);
            kx_temp = kx(:,SMS_idx_temp,iframe);
            ky_temp = ky(:,SMS_idx_temp,iframe);
            Data_cine{i}.kSpace(:,:,iframe,:,:,:,j) = GROG.GROG_Dictionary_interp_new(kSpace_temp_frame,Data{i}.G{j}.Gx,Data{i}.G{j}.Gy,kx_temp,ky_temp);
        end
    end
    Data_cine{i}.sens_map = Data{i}.sens_map;
    Data_cine{i} = get_Data_SMS(Data_cine{i},para);    
end

% Reconstrcution 
Data_cine_all.kSpace = cat(6,Data_cine{1}.kSpace,Data_cine{2}.kSpace,Data_cine{3}.kSpace);
Data_cine_all.mask = cat(6,Data_cine{1}.mask,Data_cine{2}.mask,Data_cine{3}.mask);
Data_cine_all.first_est = cat(6,Data_cine{1}.first_est,Data_cine{2}.first_est,Data_cine{3}.first_est);
Data_cine_all.sens_map = cat(6,Data_cine{1}.sens_map,Data_cine{2}.sens_map,Data_cine{3}.sens_map);
Data_cine_all.filter = Data_cine{1}.filter;

scale = max(abs(Data_cine_all.first_est(:)));
para.Recon.weight_tTV = scale*0.05;
para.Recon.weight_sTV = scale*0.001;
para.Recon.noi = 150;
Image_cine = squeeze(STCR_conjugate_gradient_MSMS(Data_cine_all,para));
%}
%% Recon systolic

% Image_sys = Image(:,:,sys,:,:);
% Data_perfusion.kSpace = Data_all.kSpace(:,:,sys,:,:,:,:);
% Data_perfusion.mask = Data_all.mask(:,:,sys,:,:,:,:);
% Data_perfusion.first_est = permute(Image(:,:,sys,:,:),[1,2,3,6,4,5]);
% Data_perfusion.sens_map = Data_all.sens_map;
% 
% scale = max(abs(Data_perfusion.first_est(:)));
% para.Recon.weight_tTV = scale*0.005;
% para.Recon.weight_sTV = scale*0.001;
% para.Recon.noi = 5;
% para.Recon.type = 'seperate SMS';
% para.Recon.no_comp = no_comp;
para.Recon.noi = 50;
sys_logical = false([1,siz(3)]);
sys_logical(sys) = true;
dia_logical = false([1,siz(3)]);
dia_logical(dia) = true;


Image_perfusion = zeros(siz,'single');

%Block_frame = 120;
%Nblock = ceil(siz(3)/Block_frame);
if Nblock>1
    Data_temp.sens_map = Data_all.sens_map;
    Data_temp.filter = Data_all.filter;
    Data_temp.SMS = Data_all.SMS;
    for i=1:Nblock
        range = (1:Block_frame) + (i-1)*Block_frame;
        if i==1
            range = [range,max(range)+1:max(range+5)];
        elseif i == Nblock
            range(range>siz(3)) = [];
            range = [min(range)-5:min(range)-1,range];
        else
            range = [min(range)-5:min(range)-1,range,max(range)+1:max(range+5)];
        end
        Data_temp.kSpace = Data_all.kSpace(:,:,range,:,:,:,:);
        Data_temp.mask = Data_all.mask(:,:,range,:,:,:,:);
        Data_temp.first_est = permute(Image(:,:,range,:,:),[1,2,3,6,4,5]);
        para.sys = sys_logical(range);
        para.dia = dia_logical(range);
        Image_temp = squeeze(STCR_conjugate_gradient_MSMS_systolic(Data_temp,para));
        if i==1
            range(end-4:end) = [];
            flag = 0;
        elseif i == Nblock
            range(1:5) = [];
            flag = 1;
        else
            range(1:5) = [];
            flag = 1;
            range(end-4:end) = [];
        end
        Image_perfusion(:,:,range,:,:) = Image_temp(:,:,range - (i-1)*Block_frame + flag*5,:,:);
    end
else
    Image_perfusion = squeeze(STCR_conjugate_gradient_MSMS_systolic(Data_all,para));
end
% save('Image_perfusion.mat','Image_perfusion','para');
% 
% range = 1:100;
% Data_temp = Data_all;
% Data_temp.kSpace = Data_all.kSpace(:,:,range,:,:,:,:);
% Data_temp.mask = Data_all.mask(:,:,range,:,:,:,:);
% Data_temp.first_est = permute(Image(:,:,range,:,:),[1,2,3,6,4,5]);
% para.sys = sys_logical(range);
% para.dia = dia_logical(range);
% Image_temp = squeeze(STCR_conjugate_gradient_MSMS_systolic(Data_temp,para));
% Image_perfusion(:,:,1:80,:,:) = Image_temp(:,:,1:80,:,:);
% 
% range = 66:165;
% Data_temp = Data_all;
% Data_temp.kSpace = Data_all.kSpace(:,:,range,:,:,:,:);
% Data_temp.mask = Data_all.mask(:,:,range,:,:,:,:);
% Data_temp.first_est = permute(Image(:,:,range,:,:),[1,2,3,6,4,5]);
% para.sys = sys_logical(range);
% para.dia = dia_logical(range);
% Image_temp = squeeze(STCR_conjugate_gradient_MSMS_systolic(Data_temp,para));
% Image_perfusion(:,:,81:160,:,:) = Image_temp(:,:,16:95,:,:);
% 
% range = 131:230;
% Data_temp = Data_all;
% Data_temp.kSpace = Data_all.kSpace(:,:,range,:,:,:,:);
% Data_temp.mask = Data_all.mask(:,:,range,:,:,:,:);
% Data_temp.first_est = permute(Image(:,:,range,:,:),[1,2,3,6,4,5]);
% para.sys = sys_logical(range);
% para.dia = dia_logical(range);
% Image_temp = squeeze(STCR_conjugate_gradient_MSMS_systolic(Data_temp,para));
% Image_perfusion(:,:,161:230,:,:) = Image_temp(:,:,31:100,:,:);

Image_sys = Image_perfusion(:,:,sys,:,:);
Image_dia = Image_perfusion(:,:,dia,:,:);

para.kSpace_info.TimeStamp_sys = para.kSpace_info.TimeStamp(sys);
para.kSpace_info.TimeStamp_dia = para.kSpace_info.TimeStamp(dia);

disp('Saving image into Results...')
if para.Recon.crop_half_FOV == 1
    Image_sys = abs(crop_half_FOV(Image_sys));
    Image_dia = abs(crop_half_FOV(Image_dia));
else
    Image_sys = abs(Image_sys);
    Image_dia = abs(Image_dia);
end

save([para.dir.save_recon_img_mat_dir,para.dir.save_recon_img_name],'Image_sys','Image_dia','Image','para','-v7.3');
disp('Reconstruction done');fprintf('\n')

% for i=1:nset
%     Data{i}.kSpace = Data{i}.kSpace(:,:,sys,:,:,:,:);
%     Data{i}.mask = Data{i}.mask(:,:,sys,:,:,:,:);
%     Data{i}.first_est = permute(Image(:,:,sys,:,i),[1,2,3,5,4]); 
% end
% 
% for i=1:nset
%     
% end

