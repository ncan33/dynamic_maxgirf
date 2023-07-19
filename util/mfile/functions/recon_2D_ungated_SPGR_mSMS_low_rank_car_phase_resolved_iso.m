function recon_2D_ungated_SPGR_mSMS_low_rank_car_phase_resolved_iso(kSpace_all,kSpace_info,para)
keyboard
%% recon PD
%kSpace_all = yExportKspace(kSpace_all,0,0.5);
set = kSpace_info.set==0;
set = find(set); set = set(1:600);
nset = max(kSpace_info.set(:))+1;
% RING trajectory correction using PD rays
para.trajectory_correction = RING_SMS(kSpace_all(:,set,:),kSpace_info.angle_mod(set));
keyboard
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
non_steady_state_rays = 360;
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
    
    if nset == 3
        nor_sl = para.Recon.nSMS*4;
    elseif nset == 4
        nor_sl = para.Recon.nSMS*3;
    elseif nset == 5
        nor_sl = 12;
    elseif nset == 6
        nor_sl = 9;
    end
    nor_sl = 12;
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

[para.kSpace_info.TimeStamp] = para.kSpace_info.TimeStamp(1:nset*nor_sl:end);

siz = size(Data{1}.first_est);
siz(4) = nset;
Image = zeros(siz,'single');

keyboard

%% recon reference image
para.setting.ifGPU = 1;
para.setting.ifplot = 0;
para.Recon.type = 'seperate SMS test';
para.Recon.noi = 50;

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
    save(fullfile(para.dir.save_recon_img_mat_dir,para.dir.save_recon_img_name),'Image','-append')
end

keyboard
% binning
% cardiac_signal = compare_curve_same_image(crop_half_FOV(Image(:,:,:,2,1)));
% cardiac_signal = normalize_to_0_1_3(smooth(cardiac_signal));
% cardiac_bins = binning(cardiac_signal,6);
% 
% respiration_signal = compare_curve_same_image(crop_half_FOV(Image(:,:,:,3,1)));
% respiration_signal = normalize_to_0_1_3(smooth(respiration_signal,30));
% respiration_bins = binning(respiration_signal,5);
% 
% para.Recon.cardiac_bins = cardiac_bins;
% para.Recon.respiration_bins = respiration_bins;
% para.Recon.bins = get_bins(cardiac_bins,respiration_bins);
%

%% self-gating and reshape k-space to resolve cardiac motion
%cardiac_signal = compare_curve_same_image(crop_half_FOV(Image(:,:,:,2,1)));
% cardaic signal
cardiac_signal = abs(crop_half_FOV(Image));
cardiac_signal = cardiac_signal(:,:,:,:);
order = vec((1:nset)'+([1,para.Recon.nSMS:-1:2]-1)*nset);
cardiac_signal = cardiac_signal(:,:,:,order);
cardiac_signal_std = std(cardiac_signal(:,:,1:round(size(Image,3)/2),:,:),1,3);
threshold = 0.5;
cardiac_signal_mask = cardiac_signal_std > max(cardiac_signal_std(:))*threshold;
while sum(cardiac_signal_mask(:)) < 0.05*numel(cardiac_signal_mask)
    threshold = threshold - 0.01;
    cardiac_signal_mask = cardiac_signal_std > max(cardiac_signal_std(:))*threshold;
end
% respiratory_signal_mask = cardiac_signal_mask;
% respiratory_signal = cardiac_signal;
cc = bwconncomp(squeeze(cardiac_signal_mask));
for i=1:length(cc.PixelIdxList)
    Npixel(i) = length(cc.PixelIdxList{i});
end
[Npixel_sort,maxN] = sort(Npixel,'descend');
Narea = sum(Npixel_sort > sum(cardiac_signal_mask(:))*0.01);

for i=1:Narea
    cardiac_signal_mask_temp = false(size(cardiac_signal_mask));
    cardiac_signal_mask_temp(cc.PixelIdxList{maxN(i)}) = true;
    center_temp = regionprops(squeeze(cardiac_signal_mask_temp),'centroid');
    centerMask(:,i) = center_temp.Centroid;
    cardiac_signal_all(:,i) = sum(sum(sum(cardiac_signal.*cardiac_signal_mask_temp,4)))/sum(cardiac_signal_mask_temp(:));
end
Distance2Center = sos(centerMask - [sx/2;sx/2;para.Recon.nSMS*nset+1]/2,1);
idx_keep = Distance2Center < sx/8;
cardiac_signal = cardiac_signal_all(:,idx_keep);

[~,cardiac_signal_max] = max(cardiac_signal);
cardiac_signal_max(cardiac_signal_max>0.6*nof) = 0;
[~,cardiac_signal_idx] = max(cardiac_signal_max);
cardiac_signal = cardiac_signal_all(:,cardiac_signal_idx);
sys = local_max(-cardiac_signal);
para.Recon.self_gating.cardiac_signal = cardiac_signal;
para.Recon.self_gating.sys = sys;

% cardiac
% 
% cardiac_signal_mask = false(size(cardiac_signal_mask));
% cardiac_signal_mask(cc.PixelIdxList{max3(2)}) = true;
% cardiac_signal_2 = sum(sum(sum(cardiac_signal.*cardiac_signal_mask,4)));
% cardiac_signal_mask = false(size(cardiac_signal_mask));
% cardiac_signal_mask(cc.PixelIdxList{max3(3)}) = true;
% cardiac_signal_3 = sum(sum(sum(cardiac_signal.*cardiac_signal_mask,4)));
% [~,loc1] = max(cardiac_signal_1);
% [~,loc2] = max(cardiac_signal_2);
% if loc1<loc2
%     cardiac_signal = cardiac_signal_2;
% else
%     cardiac_signal = cardiac_signal_1;
% end
% cardiac_signal = squeeze(cardiac_signal);


% respiratory signal
respiratory_signal = cardiac_signal_all(:,~idx_keep);
[~,respiratory_signal_idx] = max(sum(respiratory_signal));
respiratory_signal = respiratory_signal(:,respiratory_signal_idx);

% respiratory_signal_mask(cc.PixelIdxList{max2(1)}) = false;
% respiratory_signal_mask(cc.PixelIdxList{max2(2)}) = false;
% respiratory_signal = respiratory_signal.*respiratory_signal_mask;
% respiratory_signal = squeeze(sum(sum(sum(respiratory_signal)),4));
respiratory_signal = smooth(respiratory_signal,20);
para.Recon.self_gating.respiratory_signal = respiratory_signal;


% reshape k-space
% drop first and last few cardiac cycles
Image(:,:,[1:sys(2)-1,sys(end-1):end],:,:) = [];
cardiac_signal([1:sys(2)-1,sys(end-1):end]) = [];
respiratory_signal([1:sys(2)-1,sys(end-1):end]) = [];
para.kSpace_info.TimeStamp([1:sys(2)-1,sys(end-1):end]) = [];

ray_drop = [1:(sys(2)-1)*nor_sl*nset,(sys(end-1)-1)*nor_sl*nset+1:nor_all];
kSpace_all(:,ray_drop,:) = [];
theta = kSpace_info.angle_mod; theta(ray_drop) = [];
phase = kSpace_info.phase_mod; phase(ray_drop) = [];
set = kSpace_info.set+1; set(ray_drop) = [];

sys(1) = []; 
sys = sys - sys(1) + 1;
sys(end-1:end) = [];
para.kSpace_info.TimeStamp = para.kSpace_info.TimeStamp(sys);
Ncycle = length(sys); 
Nphase = median(diff(sys));

% pre-allocate arrays 
[sx,sy,nof,nset,nSMS] = size(Image);
im_bin = zeros(sx,sy,Nphase,Ncycle,nset,nSMS,'like',Image);
kSpace_bin = zeros(sx,nSMS*nset,Nphase,Ncycle,no_comp,'like',kSpace_all);
theta_bin = zeros(1,nSMS*nset,Nphase,Ncycle);
phase_bin = zeros(1,nSMS*nset,Nphase,Ncycle);
set_bin = zeros(1,nSMS*nset,Nphase,Ncycle);

dia_loc = round((Nphase+1)/2);
Ns2d = dia_loc;
Nd2s = Nphase - dia_loc + 1;
for i=1:Ncycle
    if i == Ncycle
        Frames = sys(end):nof;
        Rays = (sys(end)-1)*nor_sl*nset+1:nof*nor_sl*nset;
    else
        Frames = sys(i):sys(i+1)-1;
        Rays = (sys(i)-1)*nor_sl*nset+1:(sys(i+1)-1)*nor_sl*nset;
    end
    Image_temp = Image(:,:,Frames,:,:);
    cardiac_signal_temp = cardiac_signal(Frames);
    
    kSpace_temp = kSpace_all(:,Rays,:);
    theta_temp = theta(Rays);
    phase_temp = phase(Rays);
    set_temp = set(Rays);
    
    [~,dia_loc_current] = max(cardiac_signal_temp);
    
    im_sys2dia = Image_temp(:,:,1:dia_loc_current,:,:);
    im_dia2sys = Image_temp(:,:,dia_loc_current:end,:,:);
    
    ns2d = size(im_sys2dia,3);
    nd2s = size(im_dia2sys,3);
    
    dia_loc_rays_start = (dia_loc_current-1)*nor_sl*nset+1;
    dia_loc_rays_end = (dia_loc_current)*nor_sl*nset;
    nor_temp = length(Rays);
    
    if ns2d ~= Ns2d
        im_sys2dia = permute(im_sys2dia,[1,2,4,5,3]);
        siz = size(im_sys2dia);
        im_sys2dia = reshape(im_sys2dia,[prod(siz(1:end-1)),siz(end)]);
        [y,x] = meshgrid(1:ns2d,1:prod(siz(1:end-1)));
        [yv,xv] = meshgrid(1:(ns2d-1)/(Ns2d-1):ns2d,1:prod(siz(1:end-1)));
        im_sys2dia = interp2(y,x,im_sys2dia,yv,xv);
        im_sys2dia = reshape(im_sys2dia,[siz(1:end-1),Ns2d]);
        im_sys2dia = permute(im_sys2dia,[1,2,5,3,4]);
        im_bin(:,:,1:Ns2d,i,:,:) = im_sys2dia;
        
        rays_begin = round(1:(dia_loc_rays_start-1)/(Ns2d-1):dia_loc_rays_start);
        if ns2d<Ns2d            
            for Phase = 1:Ns2d
                kSpace_bin(:,1:nor_sl*nset,Phase,i,:) = kSpace_temp(:,rays_begin(Phase):rays_begin(Phase)+nor_sl*nset-1,:);
                theta_bin(1,1:nor_sl*nset,Phase,i) = theta_temp(1,rays_begin(Phase):rays_begin(Phase)+nor_sl*nset-1);
                phase_bin(1,1:nor_sl*nset,Phase,i) = phase_temp(1,rays_begin(Phase):rays_begin(Phase)+nor_sl*nset-1);
                set_bin(1,1:nor_sl*nset,Phase,i) = set_temp(1,rays_begin(Phase):rays_begin(Phase)+nor_sl*nset-1);
            end
        else
            nor_bin = ceil(mean(diff(rays_begin))/nSMS/nset)*nSMS*nset;
            for Phase = 1:Ns2d
                if Phase == Ns2d
                    nor_bin = nor_sl*nset;
                end
                kSpace_bin(:,1:nor_bin,Phase,i,:) = kSpace_temp(:,rays_begin(Phase):rays_begin(Phase)+nor_bin-1,:);
                theta_bin(1,1:nor_bin,Phase,i) = theta_temp(1,rays_begin(Phase):rays_begin(Phase)+nor_bin-1);
                phase_bin(1,1:nor_bin,Phase,i) = phase_temp(1,rays_begin(Phase):rays_begin(Phase)+nor_bin-1);
                set_bin(1,1:nor_bin,Phase,i) = set_temp(1,rays_begin(Phase):rays_begin(Phase)+nor_bin-1);
            end
        end
    else
        im_bin(:,:,1:Ns2d,i,:,:) = im_sys2dia;
        
        kSpace_bin(:,1:nor_sl*nset,1:Ns2d,i,:) = reshape(kSpace_temp(:,1:dia_loc_rays_end,:),[sx,nor_sl*nset,Ns2d,1,no_comp]);
        theta_bin(1,1:nor_sl*nset,1:Ns2d,i) = reshape(theta_temp(1,1:dia_loc_rays_end),[1,nor_sl*nset,Ns2d]);
        phase_bin(1,1:nor_sl*nset,1:Ns2d,i) = reshape(phase_temp(1,1:dia_loc_rays_end),[1,nor_sl*nset,Ns2d]);
        set_bin(1,1:nor_sl*nset,1:Ns2d,i) = reshape(set_temp(1,1:dia_loc_rays_end),[1,nor_sl*nset,Ns2d]);
    end
    if nd2s ~= Nd2s
        im_dia2sys = permute(im_dia2sys,[1,2,4,5,3]);
        siz = size(im_dia2sys);
        im_dia2sys = reshape(im_dia2sys,[prod(siz(1:end-1)),siz(end)]);
        [y,x] = meshgrid(1:nd2s,1:prod(siz(1:end-1)));
        [yv,xv] = meshgrid(1:(nd2s-1)/(Nd2s-1):nd2s,1:prod(siz(1:end-1)));
        im_dia2sys = interp2(y,x,im_dia2sys,yv,xv);
        im_dia2sys = reshape(im_dia2sys,[siz(1:end-1),Nd2s]);
        im_dia2sys = permute(im_dia2sys,[1,2,5,3,4]);
        im_bin(:,:,Ns2d:end,i,:,:) = im_dia2sys;
        
        rays_begin = round(dia_loc_rays_start:(nor_temp-nor_sl*nset-dia_loc_rays_start+1)/(Nd2s-1):(nor_temp-nor_sl*nset+1));
        if nd2s<Nd2s
            for Phase = 1:Nd2s
                kSpace_bin(:,1:nor_sl*nset,Phase+Ns2d-1,i,:) = kSpace_temp(:,rays_begin(Phase):rays_begin(Phase)+nor_sl*nset-1,:);
                theta_bin(1,1:nor_sl*nset,Phase+Ns2d-1,i) = theta_temp(1,rays_begin(Phase):rays_begin(Phase)+nor_sl*nset-1);
                phase_bin(1,1:nor_sl*nset,Phase+Ns2d-1,i) = phase_temp(1,rays_begin(Phase):rays_begin(Phase)+nor_sl*nset-1);
                set_bin(1,1:nor_sl*nset,Phase+Ns2d-1,i) = set_temp(1,rays_begin(Phase):rays_begin(Phase)+nor_sl*nset-1);
            end
        else
            nor_bin = ceil(mean(diff(rays_begin))/nSMS/nset)*nSMS*nset;
            for Phase = 1:Nd2s
                if Phase == Nd2s
                    nor_bin = nor_sl*nset;
                end
                kSpace_bin(:,1:nor_bin,Phase+Ns2d-1,i,:) = kSpace_temp(:,rays_begin(Phase):rays_begin(Phase)+nor_bin-1,:);
                theta_bin(1,1:nor_bin,Phase+Ns2d-1,i) = theta_temp(1,rays_begin(Phase):rays_begin(Phase)+nor_bin-1);
                phase_bin(1,1:nor_bin,Phase+Ns2d-1,i) = phase_temp(1,rays_begin(Phase):rays_begin(Phase)+nor_bin-1);
                set_bin(1,1:nor_bin,Phase+Ns2d-1,i) = set_temp(1,rays_begin(Phase):rays_begin(Phase)+nor_bin-1);
            end
        end
    else
        im_bin(:,:,Ns2d:end,i,:,:) = im_dia2sys;
        
        kSpace_bin(:,1:nor_sl*nset,Ns2d:end,i,:) = reshape(kSpace_temp(:,dia_loc_rays_start:end,:),[sx,nor_sl*nset,Nd2s,1,no_comp]);
        theta_bin(1,1:nor_sl*nset,Ns2d:end,i) = reshape(theta_temp(1,dia_loc_rays_start:end),[1,nor_sl*nset,Nd2s]);
        phase_bin(1,1:nor_sl*nset,Ns2d:end,i) = reshape(phase_temp(1,dia_loc_rays_start:end),[1,nor_sl*nset,Nd2s]);
        set_bin(1,1:nor_sl*nset,Ns2d:end,i) = reshape(set_temp(1,dia_loc_rays_start:end),[1,nor_sl*nset,Nd2s]);
    end

end

[sx,nor_bin,Nphase,Ncycle,no_comp] = size(kSpace_bin);
kSpace_bin = reshape(kSpace_bin,[sx,nor_bin,Nphase*Ncycle,no_comp]);
theta_bin = reshape(theta_bin,[1,nor_bin,Nphase*Ncycle]);
phase_bin = reshape(phase_bin,[1,nor_bin,Nphase*Ncycle]);
set_bin = reshape(set_bin,[1,nor_bin,Nphase*Ncycle]);

nof = Nphase*Ncycle;
nor_sl_new = nor_bin/nset;
for i=1:nset
    kSpace_radial = zeros([sx,nor_sl_new,nof,no_comp],'like',kSpace_bin);
    theta = zeros([1,nor_sl_new,nof],'like',theta_bin);
    phase = repmat(0:nSMS-1,[1,nor_sl_new/nSMS,nof]);
    for j=1:nof
        set = set_bin(:,:,j)==i;
        nor_temp = sum(set);
        kSpace_radial(:,1:nor_temp,j,:) = kSpace_bin(:,set,j,:);
        theta(1,1:nor_temp,j) = theta_bin(:,set,j);
        phase(1,1:nor_temp,j) = phase_bin(:,set,j);
    end

    [kx,ky] = get_k_coor(sx,theta,0,kCenter);
    
    % RING trajectory correction
    correction = para.trajectory_correction.*permute([cos(theta);sin(theta)],[4,1,2,3]);
    correction = squeeze(sum(correction,2));
    kx = kx - correction(1,:,:);
    ky = ky - correction(2,:,:);
    
    Data{i}.kSpace = GROG.SMS_GROG(kSpace_radial,kx,ky,phase,para,Data{i}.G);
    [Data{i},para] = get_Data_SMS(Data{i},para);

    para.Recon.no_comp = no_comp;
end

para.Recon.bins = repmat(diag(true(1,Nphase)),[1,Ncycle]);
im_bin = reshape(im_bin,[sx,sx,nof,nset,nSMS]);

para.setting.ifGPU = 0;
para.setting.ifplot = 0;
para.Recon.noi = 50;
save(fullfile(para.dir.save_recon_img_mat_dir,para.dir.save_recon_img_name),'para','-append')


%% reshape k-space for cine
dia = local_max(cardiac_signal);
para.Recon.self_gating.dia = dia;

cardiac_signal_dia = cardiac_signal(dia);
[~,peak_enhan] = max(cardiac_signal_dia);
[~,enhan_begin] = max(abs(diff(smooth(cardiac_signal_dia))));
first_pass_end = peak_enhan*3-enhan_begin*2;
idx = respiratory_signal(dia) < median(respiratory_signal(1:dia(first_pass_end)));
idx = idx(1:first_pass_end);
cardiac_cycles_cine = find(diff(find(idx))==1);
idx = find(idx);
cardiac_cycles_cine = idx(cardiac_cycles_cine);
Ncycles_cine = length(cardiac_cycles_cine);

theta = kSpace_info.angle_mod; theta(ray_drop) = [];
phase = kSpace_info.phase_mod; phase(ray_drop) = [];
set = kSpace_info.set+1; set(ray_drop) = [];

kSpace_cine = zeros(sx,1,Ncycles_cine,no_comp,nset);
theta_cine = zeros(1,1,Ncycles_cine,nset);
phase_cine = zeros(1,1,Ncycles_cine,nset);

for i=1:Ncycles_cine
    cycle_begin_ray = (dia(cardiac_cycles_cine(i))-1)*nor_sl+1;
    cycle_end_ray = (dia(cardiac_cycles_cine(i)+1)-1)*nor_sl;
    selected_rays = cycle_begin_ray:cycle_end_ray;
    nor_temp = floor(length(selected_rays)/nSMS)*nSMS;
    selected_rays(nor_temp+1:end) = [];
    for j=1:nset
        set_temp = find(set==j);
        idx_temp = set_temp(selected_rays);
        kSpace_temp = kSpace_all(:,idx_temp,:);
        theta_temp = theta(idx_temp);
        phase_temp = phase(idx_temp);
        
        kSpace_cine(:,1:nor_temp,i,:,j) = kSpace_temp;
        theta_cine(1,1:nor_temp,i,j) = theta_temp;
        phase_cine(1,1:nor_temp,i,j) = phase_temp;
    end
end

nor_cine_each_cycle = squeeze(sum(kSpace_cine(sx/2,:,:,1,1)~=0));
nor_cine = floor(median(nor_cine_each_cycle)/nSMS)*nSMS;
cine_window_length = nSMS;
cine_sliding_length = nSMS;
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
Data_cine_all = Data_cine{1};
for i=2:nset
    Data_cine_all.kSpace = cat(6,Data_cine_all.kSpace,Data_cine{i}.kSpace);
    Data_cine_all.mask = cat(6,Data_cine_all.mask,Data_cine{i}.mask);
    Data_cine_all.first_est = cat(6,Data_cine_all.first_est,Data_cine{i}.first_est);
    Data_cine_all.sens_map = cat(6,Data_cine_all.sens_map,Data_cine{i}.sens_map);
end
for i=1:para.Recon.no_comp
    k_temp = Data_cine_all.kSpace(:,:,:,i,:,:,:);
    kSpace(:,i) = k_temp(Data_cine_all.mask);
end
Data_cine_all.kSpace = kSpace;
clear kSpace k_temp

scale = max(abs(Data_cine_all.first_est(:)));
para.Recon.weight_tTV = scale*0.05;
para.Recon.weight_sTV = scale*0.001;
para.Recon.type = 'seperate SMS test';
para.Recon.noi = 150;

clearvars -except Data Data_cine_all para im_bin nSMS nset

Image_cine = squeeze(STCR_conjugate_gradient_MSMS_cine(Data_cine_all,para));
Image_cine = abs(crop_half_FOV(Image_cine));
Image_cine = Image_cine(:,:,:,:);
order = vec([1:nSMS:nSMS*nset;(1:nSMS:nSMS*nset)+(nSMS-1:-1:1)']');
Image_cine = Image_cine(:,:,:,order);
Image_cine = permute(Image_cine,[1,2,4,3]);
save(fullfile(para.dir.save_recon_img_mat_dir,para.dir.save_recon_img_name),'Image_cine','-append')

%% 2nd reconstruction stage

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

Data_all.first_guess = permute(im_bin,[1,2,3,6,5,4]);
para.Recon.weight_tTV = scale*0.04;
para.Recon.weight_sTV = scale*0.00;
para.Recon.type = 'seperate SMS test';
para.Recon.noi = 50;

clearvars -except Data_all para

[Data_all,para] = get_iso_image(Data_all,para);
Data_all = Patch_tracking(Data_all,para);

im_bin = STCR_conjugate_gradient_MSMS_iso_patch_llr(Data_all,para);
save(fullfile(para.dir.save_recon_img_mat_dir,para.dir.save_recon_img_name),'im_bin','para','-append')

Image_sys = abs(crop_half_FOV(im_bin(:,:,para.Recon.bins(1,:),:)));
Image_sys = permute(Image_sys,[1,2,4,3]);
save(fullfile(para.dir.save_recon_img_mat_dir,para.dir.save_recon_img_name),'Image_sys','-append')

dia_loc = round((size(para.Recon.bins,1)+1)/2);
Image_dia = abs(crop_half_FOV(im_bin(:,:,para.Recon.bins(dia_loc,:),:)));
Image_dia = permute(Image_dia,[1,2,4,3]);
save(fullfile(para.dir.save_recon_img_mat_dir,para.dir.save_recon_img_name),'Image_dia','-append')

return












%{
%%
% for i=1:nset
%     scale = max(abs(Data{i}.first_est(:)));
%     para.Recon.weight_tTV = scale*0.04;
%     para.Recon.weight_sTV = scale*0.00;
%     for j=1:para.Recon.no_comp
%         k_temp = Data{i}.kSpace(:,:,:,j,:,:,:);
%         kSpace(:,j) = k_temp(Data{i}.mask);
%     end
%     Data{i}.kSpace = kSpace;
%     clear kSpace
%     Data{i}.first_guess = im_bin(:,:,:,i,:);
%     im_bin(:,:,:,i,:) = STCR_conjugate_gradient_low_rank_bins(Data{i},para);
%     save(fullfile(para.dir.save_recon_img_mat_dir,para.dir.save_recon_img_name),'im_bin','-append')
% end


% 
% 
% Data_all = Data{1};
% for i=2:nset
%     Data_all.kSpace = cat(6,Data_all.kSpace,Data{i}.kSpace);
%     Data_all.mask = cat(6,Data_all.mask,Data{i}.mask);
%     Data_all.first_est = cat(6,Data_all.first_est,Data{i}.first_est);
%     Data_all.sens_map = cat(6,Data_all.sens_map,Data{i}.sens_map);
% end
% clear Data
% 
% 
% for i=1:para.Recon.no_comp
%     k_temp = Data_all.kSpace(:,:,:,i,:,:,:);
%     kSpace(:,i) = k_temp(Data_all.mask);
% end
% Data_all.kSpace = kSpace;
% clear kSpace k_temp
% 
% scale = max(abs(Data_all.first_est(:)));
% para.Recon.weight_tTV = scale*0.04;
% para.Recon.weight_sTV = scale*0.001;
% para.Recon.type = 'seperate SMS test';
% para.Recon.noi = 50;
% 
% Image = STCR_conjugate_gradient_MSMS_low_rank(Data_all,para);
% % recon by blocks
% siz = size(squeeze(Data_all.first_est));
% Image = zeros(siz,'single');
% 
% Block_frame = 70;
% Nblock = ceil(siz(3)/Block_frame);
% 
% if Nblock>1
%     Data_temp.sens_map = Data_all.sens_map;
%     Data_temp.filter = Data_all.filter;
%     Data_temp.SMS = Data_all.SMS;
%     for i=1:Nblock
%         range = (1:Block_frame) + (i-1)*Block_frame;
%         if i==1
%             range = [range,max(range)+1:max(range+5)];
%         elseif i == Nblock
%             range(range>siz(3)) = [];
%             range = [min(range)-5:min(range)-1,range];
%         else
%             range = [min(range)-5:min(range)-1,range,max(range)+1:max(range+5)];
%         end
%         Data_temp.kSpace = Data_all.kSpace(:,:,range,:,:,:,:);
%         Data_temp.mask = Data_all.mask(:,:,range,:,:,:,:);
%         Data_temp.first_est = Data_all.first_est(:,:,range,:,:,:,:);
%         Image_temp = squeeze(STCR_conjugate_gradient_MSMS(Data_temp,para));
%         if i==1
%             range(end-4:end) = [];
%             flag = 0;
%         elseif i == Nblock
%             range(1:5) = [];
%             flag = 1;
%         else
%             range(1:5) = [];
%             flag = 1;
%             range(end-4:end) = [];
%         end
%         Image(:,:,range,:,:) = Image_temp(:,:,range - (i-1)*Block_frame + flag*5,:,:);
%     end
% else
%     Image = squeeze(STCR_conjugate_gradient_MSMS(Data_all,para));
% end
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
%}
      