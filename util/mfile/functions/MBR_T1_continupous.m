function [Image,T1_map,para] = MBR_T1_continupous(Data,para)

disp('Performing iterative STCR reconstruction...');
disp('Showing progress...')


ifplot         = para.setting.ifplot;
ifGPU          = para.setting.ifGPU;
weight_tTV     = para.Recon.weight_tTV;
weight_sTV     = para.Recon.weight_sTV;
beta_sqrd      = para.Recon.epsilon;
para.Recon.step_size_1 = para.Recon.step_size(1);

nor = para.Recon.nor;
% nof_all = para.nor_total/nor;
% nof_one = nof_all/para.nof_new;

new_img_x_1 = single(Data.first_est);


if isfield(Data,'sens_map')
    Data.sens_map_conj = conj(Data.sens_map);
end

if isfield(Data,'phase_mod')
    Data.phase_mod_conj = conj(Data.phase_mod);
end

if ifGPU
    Data.kSpace        = gpuArray(Data.kSpace);
    new_img_x_1          = gpuArray(new_img_x_1);
    Data.sens_map      = gpuArray(Data.sens_map);
    Data.sens_map_conj = gpuArray(Data.sens_map_conj);
    if isfield(Data,'mask')
        Data.mask          = gpuArray(Data.mask);
    end
    if isfield(Data,'filter')
        Data.filter        = gpuArray(Data.filter);
    end
    
    gpuInfo = gpuDevice;
    gpuSize = gpuInfo.AvailableMemory;
    imSize  = numel(new_img_x_1)*8;
    if imSize*para.Recon.no_comp > gpuSize*0.3
        para.Recon.type = [para.Recon.type,' less memory'];
    end
end

para.Cost_1 = struct('fidelityNorm',[],'temporalNorm',[],'spatialNorm',[],'totalCost',[]);

fidelity = @(im) compute_fidelity_yt_new(im,Data,para);
spatial  = @(im) compute_sTV_yt(im,weight_sTV,beta_sqrd);
temporal = @(im) compute_tTV_yt(im,weight_tTV,beta_sqrd);

%%%%% dictionary
fa = 0.4:0.05:1.2;
T1_map = 1:10:3500;
flip_angle = 6;

TI1 = 10.98;
TR = 4.76;
clear Dic
for iT1=1:length(T1_map)
    for ifa = 1:length(fa)
        Dic(:,iT1,ifa) = Bloch_IR_continupous(TI1,0.92,T1_map(iT1),[20,2600/20],flip_angle*fa(ifa),2600,TR);
    end
end
Dic = Dic(:,:);


% dweight = weight_tTV/10;
% patch_siz = [5,5];
% over_lapping = [2,2];
% for i=1:size(new_img_x_1,5)
%     P{i} = Patch_search(new_img_x_1(:,:,2,:,i),patch_siz,over_lapping);
% end
for iter_no = 1:para.Recon.noi
%%%%% fidelity term/temporal/spatial TV
%     if iter_no<=10
%         weight_tTV = weight_tTV - dweight;
%         para.Recon.weight_tTV = para.Recon.weight_tTV - dweight;
%     end

    %new_img_x_1 = LowRank_yt(new_img_x_1);

    [update_term,fidelity_norm] = compute_fidelity_yt_new(new_img_x_1,Data,para);
    [update_term_Model,T1_map] = T1_fitting_SMS_2_SLG_pattern_recognition(LowRank_yt(new_img_x_1),Dic);

%     phase = exp(1i*angle(new_img_x_1));
%     for j=1:size(phase,5)
%         new_img_x_1(:,:,:,:,j) = registration_2D_more2more(abs(new_img_x_1(:,:,:,:,j)),abs(update_term_Model(:,:,:,:,j))).*phase(:,:,:,:,j);
%     end
    update_term = update_term + (update_term_Model - new_img_x_1)*0.8;
    
%     for i=1:size(new_img_x_1,5)
%         update_term_llr(:,:,:,:,i) = Patch_low_rank_multi_contrast(new_img_x_1(:,:,:,:,i),patch_siz,over_lapping,P{i});
%     end
    
%     update_term_1 = update_term_1 + (update_term_llr - new_img_x_1)*0.5;
%     update_term = update_term + compute_tTV_yt(new_img_x_1,weight_tTV,beta_sqrd);
    update_term = update_term + spatial(new_img_x_1);

%%%%% conjugate gradient
    if iter_no > 1
        beta_1 = update_term(:)'*update_term(:)/(update_term_old_1(:)'*update_term_old_1(:)+eps('single'));
        update_term = update_term + beta_1*update_term_old_1;
    end
    update_term_old_1 = update_term; clear update_term_1

    
%%%%% line search  
    
%     para.Cost = Cost_STCR(fidelity_norm, new_img_x, weight_sTV, weight_tTV, Data,update_term_M,para.Cost); clear fidelity_update
%     step_size = line_search(new_img_x,update_term_old,Data,para,update_term_M);
%     para.Recon.step_size(iter_no) = step_size;
% 
%     new_img_x = new_img_x + step_size * update_term_old;
    
    para.Cost_1 = Cost_STCR(fidelity_norm, new_img_x_1, weight_sTV, weight_tTV, Data,update_term_Model, para.Cost_1); clear fidelity_update_1
    step_size_1 = line_search(new_img_x_1,update_term_old_1,Data,para,update_term_Model,para.Cost_1.totalCost(end));
    para.Recon.step_size_1(iter_no) = step_size_1;
    new_img_x_1 = new_img_x_1 + step_size_1 * update_term_old_1;
    

%%%%% plot&save part 

    if ifplot == 1
        showImage(new_img_x_1,para.Cost_1)
        figure(2)
        im1 = crop_half_FOV(T1_map);
        imagesc([im1(:,:)])
        axis image
        drawnow
    end
end

Image = squeeze(gather(new_img_x_1));

end


function [Cost_new,Cost,fNorm,tNorm,sNorm] = Cost_STCR(fUpdate, Image, sWeight, tWeight, Data, update_term_M, Cost_old)

%N = numel(Image);

fNorm = sum(abs(fUpdate(:)).^2)*10;
mNorm = sum(vec(abs(crop_half_FOV(Image - update_term_M))).^2);

Image = crop_half_FOV(Image);

if tWeight ~= 0
    tNorm = tWeight .* abs(diff(Image,1,3));
    tNorm = sum(tNorm(:));
else
    tNorm = 0;
end

%Image = Image - crop_half_FOV(Data.ref);

if sWeight ~= 0
    sx_norm = abs(diff(Image,1,2));
    sx_norm(:,end+1,:,:,:)=0;
    sy_norm = abs(diff(Image,1,1));
    sy_norm(end+1,:,:,:,:)=0;
    sNorm = sWeight .* sqrt(abs(sx_norm).^2+abs(sy_norm).^2);
    sNorm = sum(sNorm(:));
else
    sNorm = 0;
end

Cost = sNorm + tNorm + fNorm + mNorm;
tNorm = tNorm + mNorm;

if nargin == 6
    Cost_new = Cost;
    return
end

Cost_new = Cost_old;

if isempty(Cost_old.fidelityNorm)==1
    Cost_new.fidelityNorm = gather(fNorm);
    Cost_new.temporalNorm = gather(tNorm);
    Cost_new.spatialNorm = gather(sNorm);
    Cost_new.totalCost = gather(Cost);
else    
    Cost_new.fidelityNorm(end+1) = gather(fNorm);
    Cost_new.temporalNorm(end+1) = gather(tNorm);
    Cost_new.spatialNorm(end+1) = gather(sNorm);
    Cost_new.totalCost(end+1) = gather(Cost);
end

end


function step = line_search(old,update,Data,para,update_term_M,cost_old)

step_start = para.Recon.step_size(end)*1.3;%magic number
tau = 0.8;
max_try = 15;
step = step_start;

for i=1:max_try
    
    new = old + step*update;
    fidelity_new = compute_fidelity_for_line_search_yt(new,Data,para);
    cost_new = Cost_STCR(fidelity_new,new,para.Recon.weight_sTV,para.Recon.weight_tTV,Data,update_term_M);

    if cost_new > cost_old
        step = step*tau;
    else
        return
    end

end

end


function [Image,T1_map_all] = T1_fitting_SMS_2_SLG_pattern_recognition(Image,Dic)

nof = size(Image,3);

siz = [8*150/nof,nof/8];

fa = 0.4:0.05:1.2;
T1 = 1:10:3500;

if length(vec(Dic)) < 1000
    InversionTime = Dic;
    clear Dic
    for iT1=1:length(T1)
        for ifa = 1:length(fa)
            Dic(:,iT1,ifa) = Bloch_2_IR_1_SLG(InversionTime,0.92,T1(iT1),siz,flip_angle*fa(ifa),InversionTime(6)*1000-100);
        end
    end
    Dic = Dic(1:nof,:);
end
Dic_norm = sqrt(sum(Dic.^2));
%Dic_norm = min(Dic,[],1);
Dic = Dic./Dic_norm;
Dic = permute(Dic,[3,1,2]);

Image = squeeze(Image);
im_0 = abs(Image);
% non local mean filtering
% 
% for i=1:nof
%     for j=1:size(im_0,4)
%         im_0(:,:,i,j) = imnlmfilt(im_0(:,:,i,j),'DegreeOfSmoothing',0.1);
%     end
% end
%%
im_p = angle(Image);
im_p = exp(1i.*im_p);
im_norm = sqrt(sum(im_0.^2,3));
%im_norm = min(im_0,[],3);
im_0 = im_0./im_norm;

siz = size(im_0);
if length(siz)<4
    siz(4) = 1;
end

%% fit
% T1_map_all = abs(zeros(siz(1),siz(2),siz(4),'like',Image));

tic
for i=1:siz(4)
    for j=1:siz(2)
        im = im_0(:,j,:,i);
        im = reshape(im,[siz(1),siz(3)]);
        
        d0 = im - abs(Dic);
        %d0(:,[1],:) = 0;
        %d0(:,1) = d0(:,1)/2;
        %d0(:,6) = d0(:,6)/2;
%        d0(:,[1:3,8:10,11:13,18:20,21:23,28:30,31:33,38:40,41:43,48:50],:,:) = [];
        d0 = sqrt(squeeze(sum(d0.^2,2)));
        
        [~,idx_T1] = min(min(d0,[],3),[],2);
        
        idx_map_all(:,j,i) = idx_T1;
    end
end  
toc

%

[T1_map_all,FA_map_all] = ind2sub([length(T1),length(fa)],idx_map_all);

% spatial constraint

% im_norm = im_norm + compute_sTV_yt(im_norm,3,1e-7);
% T1_map_all = T1_map_all + compute_sTV_yt(T1_map_all,20,1e-7);

%T1_map_all = round(T1_map_all);
%T1_map_all(T1_map_all>length(T1)) = length(T1);
%T1_map_all(T1_map_all<1) = 1;

idx_map_all = sub2ind([length(T1),length(fa)],T1_map_all,FA_map_all);
MBI = Dic(:,1:nof,round(idx_map_all));

T1_map_all = T1(T1_map_all);
FA_map_all = fa(FA_map_all);


MBI = reshape(MBI,[nof,siz([1,2,4])]);
MBI = permute(MBI,[2,3,1,4]);
MBI = MBI.*im_norm;

%im_p = (mean(im_p,3)-im_p)*0.9+im_p;
%im_p = im_p + compute_tTV_yt(im_p,0.1,1e-7);
MBI = MBI(:,:,1:nof,:).*im_p;

%MBI = gpuArray(MBI);

%Image(73:216,73:216,:,:) = MBI; 
Image = permute(MBI,[1,2,3,5,4]);
end