function [Image_PD,Data] = recon_PD_multiple_frames_low_rank(kSpace_all,kSpace_info,para)

kSpace_info.angle_mod = kSpace_info.angle_mod(:,1:kSpace_info.NumberOfPDlines,:);
kSpace_info.phase_mod = kSpace_info.phase_mod(:,1:kSpace_info.NumberOfPDlines,:);
kSpace_info.set = kSpace_info.set(:,1:kSpace_info.NumberOfPDlines,:);

para.kSpace_info = kSpace_info;

[sx,nor_all,no_comp] = size(kSpace_all);

para.Recon.nSMS = max(kSpace_info.phase_mod(:))+1;
kCenter = para.kSpace_center;
para.setting.ifplot = 1;

nset = max(kSpace_info.set(:))+1;

%% recon reference 
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
    
    [Data{i}.kSpace,Data{i}.G] = GROG.SMS_GROG(kSpace_radial,kx,ky,phase,para);
    para.Recon.nor = nor_sl;
    [Data{i},para] = get_Data_SMS(Data{i},para);

    para.Recon.no_comp = no_comp;
end

siz = size(Data{1}.first_est);
siz(4) = nset;
Image_PD = zeros(siz,'single');
para.setting.ifGPU = 1;

for i=1:nset
    scale = max(abs(Data{i}.first_est(:)));
    para.Recon.weight_tTV = scale*0.05;
    para.Recon.weight_sTV = scale*0.008;
    para.setting.ifplot = 0;
    para.Recon.type = 'seperate SMS test';
    para.Recon.noi = 50;
    for j=1:para.Recon.no_comp
        k_temp = Data{i}.kSpace(:,:,:,j,:,:,:);
        kSpace(:,j) = k_temp(Data{i}.mask);
    end
    Data{i}.kSpace = kSpace;
    Image_PD(:,:,:,i,:) = STCR_conjugate_gradient_low_rank_bins(Data{i},para);
    clear kSpace
end


keyboard


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
    
    [Data{i}.kSpace,Data{i}.G] = GROG.SMS_GROG(kSpace_radial,kx,ky,phase,para);
    para.Recon.nor = nor_sl;
    [Data{i},para] = get_Data_SMS(Data{i},para);

    para.Recon.no_comp = no_comp;
end


Data_all = rmfield(Data{1},{'G','first_est'});
scale = max(abs(Data{1}.first_est(:)));
for i=2:nset
    Data_all.kSpace = cat(6,Data_all.kSpace,Data{i}.kSpace);
    Data_all.mask = cat(6,Data_all.mask,Data{i}.mask);
    scale = max([scale;abs(Data{i}.first_est(:))]);
    Data_all.sens_map = cat(6,Data_all.sens_map,Data{i}.sens_map);
end
clear Data
for i=1:para.Recon.no_comp
    k_temp = Data_all.kSpace(:,:,:,i,:,:,:);
    kSpace(:,i) = k_temp(Data_all.mask);
end
Data_all.kSpace = kSpace;
clear kSpace k_temp

Data_all.first_guess = permute(Image_PD,[1,2,3,6,5,4]);
para.Recon.weight_tTV = scale*0.04;
para.Recon.weight_sTV = scale*0.00;
para.Recon.type = 'seperate SMS test';
para.Recon.noi = 50;

% Data_all.myo_mask = mask;
clearvars -except Data_all para


[sx,sy,nof,~,nSMS,nset] = size(Data_all.first_guess);
nslice = nSMS*nset;
order = 1:nSMS:nslice;
for i=nSMS:-1:2
    order = [order,i:nSMS:nslice];
end      
[~,order_back] = sort(order);

Data_all.first_guess = reshape(Data_all.first_guess,[sx,sy,nof,nslice]);
Data_all.first_guess = Data_all.first_guess(:,:,:,order);



Data_all.iso.Nphase = 1;
Data_all.iso.Ncycle = nof/Data_all.iso.Nphase;
Data_all.iso.order = order;
Data_all.iso.order_back = order_back;
Data_all.iso.nSMS = nSMS;
Data_all.iso.nset = nset;


clearvars -except Data_all para

patch_size = [5,5,2];
search_size = [9,9,4];
patch_shift = [5,5,2];
Nphase = 1;
Data_all.shifts = zeros(3,size(Data_all.first_guess,3));
Data_all.llr = Patch_tracking_3D_with_guide_random_search(Data_all.first_guess,Nphase,Data_all.shifts,patch_size,search_size,patch_shift);


para.setting.ifGPU = 1;
para.Recon.noi = 50;

para.Recon.bins = true(1,50);

im_bin = LLR_patch_tracking_4D_ADMM(Data_all,para);





Image_PD = abs(crop_half_FOV(Image_PD(:,:,:,:)));
Image_PD = squeeze(Image_PD);
order = vec((1:nset)'+([1,para.Recon.nSMS:-1:2]-1)*nset);
Image_PD = Image_PD(:,:,:,order);
