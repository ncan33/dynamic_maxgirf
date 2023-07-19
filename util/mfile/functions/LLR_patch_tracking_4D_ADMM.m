function Image = LLR_patch_tracking_4D_ADMM(Data,para)

disp('Performing ADMM reconstruction...');
disp('Showing progress...')

para.setting.ifGPU = 1;
para.Recon.noi = 40;


% Data.llr = llr_no_overlapping;
Image = STCR_conjugate_gradient_MSMS_ADMM1(Data,para);

Data.Y = zeros(size(Image));

lr_weight = 0.01;

for i=1:2
    
    Image_mc_lr = patch_low_rank_only(Image,Data.llr);

    Image_lr = reshape(Image,[para.Recon.sx,para.Recon.sx,Data.iso.Nphase,Data.iso.Ncycle,Data.iso.nset*Data.iso.nSMS]);
    Image_lr = permute(Image_lr,[1,2,5,4,3]);
    Image_lr = LowRank_patch_yt(Image_lr);
    Image_lr = permute(Image_lr,[1,2,5,4,3]);
    Image_lr = reshape(Image_lr,[para.Recon.sx,para.Recon.sx,Data.iso.Nphase*Data.iso.Ncycle,Data.iso.nset*Data.iso.nSMS]);
    Image_lr = permute(Image_lr,[1,2,4,3]);
    
    Image_mc_lr = permute(Image_mc_lr,[1,2,4,3]);
    Image_lr(Data.llr.mask) = Image_mc_lr(Data.llr.mask);
    
    para.Recon.noi = 30;
    Image_lr = permute(Image_lr,[1,2,4,3]);
    Data.first_guess = Image_lr;
    
    Data.Y = Data.Y + Image - Image_lr;

    Image = STCR_conjugate_gradient_MSMS_ADMM2(Data,para,lr_weight);
end


end


function [Image,para] = STCR_conjugate_gradient_MSMS_ADMM1(Data,para)
%[Image,para] = STCR_conjugate_gradient(Data,para)
disp('Performing iterative STCR reconstruction...');
disp('Showing progress...')

ifplot         = para.setting.ifplot;
ifGPU          = para.setting.ifGPU;
weight_tTV     = para.Recon.weight_tTV;  
weight_sTV     = para.Recon.weight_sTV;
beta_sqrd      = para.Recon.epsilon;
para.Recon.step_size = para.Recon.step_size(1);

if isfield(Data,'first_guess')
    new_img_x = Data.first_guess;   
else
    new_img_x = single(Data.first_est);
end

if isfield(Data,'phase_mod')
    Data.phase_mod_conj = conj(single(Data.phase_mod));
end

if isfield(Data,'sens_map')
    Data.sens_map_conj = conj(Data.sens_map);
end

if ifGPU
%    Data.kSpace        = gpuArray(Data.kSpace);
    new_img_x          = gpuArray(new_img_x);
    Data.sens_map      = gpuArray(Data.sens_map);
    Data.sens_map_conj = gpuArray(Data.sens_map_conj);
%     if isfield(Data,'mask')
%         Data.mask          = gpuArray(Data.mask);
%     end
    if isfield(Data,'filter')
        Data.filter        = gpuArray(Data.filter);
    end
    %Data.first_est = gpuArray(Data.first_est);
%    Data.phase_mod = gpuArray(Data.phase_mod);
%    Data.phase_mod_conj = gpuArray(Data.phase_mod_conj);
%    beta_sqrd = gpuArray(beta_sqrd);
    if isfield(Data,'N')
        for i=1:length(Data.N)
            Data.N(i).S = gpuArray(Data.N(i).S);
            Data.N(i).Apodizer = gpuArray(Data.N(i).Apodizer);
            Data.N(i).W = gpuArray(Data.N(i).W);
        end
    end
    
%     gpuInfo = gpuDevice;
%     gpuSize = gpuInfo.AvailableMemory;
%     imSize  = numel(new_img_x)*8;
%     if imSize*para.Recon.no_comp > gpuSize*0.3
%        para.Recon.type = [para.Recon.type,' less memory'];
%     end
end

para.Cost = struct('fidelityNorm',[],'temporalNorm',[],'spatialNorm',[],'totalCost',[]);
%new_img_x = new_img_x./Data.map;

fidelity = @(im) compute_fidelity_yt_new(im,Data,para);
% spatial  = @(im) compute_sTV_3D_iso(im,weight_sTV,beta_sqrd);
temporal = @(im) compute_tTV_yt(im,weight_tTV,beta_sqrd);

for iter_no = 1:para.Recon.noi

%     if mod(iter_no,1) == 1
        t1 = tic;
%     end

%%%%% fidelity term/temporal/spatial TV

    [update_term,fidelity_norm] = fidelity(iso_forward(new_img_x,Data.iso));
    update_term = iso_backward(update_term,Data.iso);
    update_term = update_term + temporal(new_img_x)*0.5;
    
    update_patch = patch_ttv(gather(new_img_x),Data.llr,para)*0.5;
    update_bins = compute_tTV_bins_no_lr(gather(new_img_x),weight_tTV,beta_sqrd,para.Recon.bins)*0.5;
    update_bins(permute(Data.llr.mask,[1,2,4,3])) = update_patch(permute(Data.llr.mask,[1,2,4,3]));
    update_term = update_term + update_bins; clear update_bins update_patch
%     update_term = update_term + tv_bin(new_img_x,Data,para)*0.5;

%%%%% conjugate gradient
    tic;
    if iter_no > 1
        beta = update_term(:)'*update_term(:)/(update_term_old(:)'*update_term_old(:)+eps('single'));
        update_term = update_term + beta*update_term_old;
    end
    update_term_old = update_term; clear update_term
    
%%%%% line search    
    
    para.Cost = Cost_STCR(fidelity_norm, new_img_x, weight_sTV, weight_tTV, para.Cost); clear fidelity_update
    step_size = line_search_iso(new_img_x,update_term_old,Data,para);
    para.Recon.step_size(iter_no) = step_size;

    new_img_x = new_img_x + step_size * update_term_old;
    para.CPUtime.update(iter_no) = toc;

%%%%% plot&save part 

    if ifplot == 1
        showImage3D(new_img_x,para.Cost)
    end

%%%%% break when step size too small or cost not changing too much

    if para.Recon.break && iter_no > 1
        if step_size<1e-4 %|| abs(para.Cost.totalCost(end) - para.Cost.totalCost(end-1))/para.Cost.totalCost(end-1) < 1e-4
            break
        end
    end
    
%     if mod(iter_no,1) == 0
        fprintf(['Iteration = ' num2str(iter_no) '...']);
        toc(t1);
%     end
end

Image = squeeze(gather(new_img_x));
%para = get_CPU_time(para);
%fprintf(['Iterative STCR running time is ' num2str(para.CPUtime.interative_recon) 's' '\n'])
end




function Image = iso_forward(Image,iso)
[sx,sy,nof,Nslice] = size(Image);
Image = Image(:,:,:,iso.order_back);
Image = reshape(Image,[sx,sy,nof,1,iso.nSMS,iso.nset]);
end

function Image = iso_backward(Image,iso)

[sx,sy,nof,~,~,~] = size(Image);
Image = reshape(Image,[sx*sy*nof,iso.nSMS*iso.nset]);
Image = Image(:,iso.order);
Image = reshape(Image,[sx,sy,nof,iso.nSMS*iso.nset]);

end


function step = line_search_iso(old,update,Data,para)

step_start = para.Recon.step_size(end)*1.3;%magic number
%step_start = 2;
%step_start = para.Recon.step_size(1);
tau = 0.8;
max_try = 15;
step = step_start;

cost_old = para.Cost.totalCost(end);

for i=1:max_try
    
    new = old + step*update;
    fidelity_new = compute_fidelity_for_line_search_yt(iso_forward(new,Data.iso),Data,para);
    switch para.Recon.type
        case '3D'
            cost_new = Cost_STCR_3D(fidelity_new,new,para.Recon.weight_sTV,para.Recon.weight_tTV,para.Recon.weight_sliceTV);
        otherwise
            cost_new = Cost_STCR(fidelity_new,new,para.Recon.weight_sTV,para.Recon.weight_tTV);
    end
    
    if cost_new > cost_old
        step = step*tau;
    else
        %fprintf(['Step = ' num2str(step) '...\n'])
        %fprintf(['Cost = ' num2str(round(cost_new)) '...\n'])
        return
    end

end
%fprintf(['Step = ' num2str(step) '...\n'])
%fprintf(['Cost = ' num2str(round(cost_new)) '...\n'])
end


% function update_tv_all = tv_bin(Image,Data,para)
% 
% [sx,sy,nof,nslice] = size(Image);
% Nphase = Data.iso.Nphase;
% Ncycle = Data.iso.Ncycle;
% update_tv_all = compute_3DtTV_yt(reshape(Image,[sx,sy,Nphase,Ncycle,nslice]),para.Recon.weight_tTV,1e-7);
% update_tv_all = reshape(update_tv_all,[sx,sy,nof,nslice]);
% 
% end



function [Image,para] = STCR_conjugate_gradient_MSMS_ADMM2(Data,para,lr_weight)
%[Image,para] = STCR_conjugate_gradient(Data,para)
disp('Performing iterative STCR reconstruction...');
disp('Showing progress...')

ifplot         = para.setting.ifplot;
ifGPU          = para.setting.ifGPU;
weight_tTV     = para.Recon.weight_tTV;
weight_sTV     = para.Recon.weight_sTV;
beta_sqrd      = para.Recon.epsilon;
para.Recon.step_size = para.Recon.step_size(1);

if isfield(Data,'first_guess')
    new_img_x = Data.first_guess;   
else
    new_img_x = single(Data.first_est);
end

if isfield(Data,'phase_mod')
    Data.phase_mod_conj = conj(single(Data.phase_mod));
end

if isfield(Data,'sens_map')
    Data.sens_map_conj = conj(Data.sens_map);
end

if ifGPU
%    Data.kSpace        = gpuArray(Data.kSpace);
    new_img_x          = gpuArray(new_img_x);
    Data.sens_map      = gpuArray(Data.sens_map);
    Data.sens_map_conj = gpuArray(Data.sens_map_conj);
%     if isfield(Data,'mask')
%         Data.mask          = gpuArray(Data.mask);
%     end
    if isfield(Data,'filter')
        Data.filter        = gpuArray(Data.filter);
    end
    %Data.first_est = gpuArray(Data.first_est);
%    Data.phase_mod = gpuArray(Data.phase_mod);
%    Data.phase_mod_conj = gpuArray(Data.phase_mod_conj);
%    beta_sqrd = gpuArray(beta_sqrd);
    if isfield(Data,'N')
        for i=1:length(Data.N)
            Data.N(i).S = gpuArray(Data.N(i).S);
            Data.N(i).Apodizer = gpuArray(Data.N(i).Apodizer);
            Data.N(i).W = gpuArray(Data.N(i).W);
        end
    end
    
%     gpuInfo = gpuDevice;
%     gpuSize = gpuInfo.AvailableMemory;
%     imSize  = numel(new_img_x)*8;
%     if imSize*para.Recon.no_comp > gpuSize*0.3
%        para.Recon.type = [para.Recon.type,' less memory'];
%     end
end

para.Cost = struct('fidelityNorm',[],'temporalNorm',[],'spatialNorm',[],'totalCost',[]);
%new_img_x = new_img_x./Data.map;

fidelity = @(im) compute_fidelity_yt_new(im,Data,para);
% spatial  = @(im) compute_sTV_3D_iso(im,weight_sTV,beta_sqrd);
% temporal = @(im) compute_tTV_yt(im,weight_tTV,beta_sqrd);

for iter_no = 1:para.Recon.noi

%     if mod(iter_no,1) == 1
        t1 = tic;
%     end

%%%%% fidelity term/temporal/spatial TV

    [update_term,fidelity_norm] = fidelity(iso_forward(new_img_x,Data.iso));
    update_term = iso_backward(update_term,Data.iso);
%    update_term = update_term + temporal(new_img_x)*0.5;
%     update_term = update_term + tv_bin(new_img_x,Data,para)*0.5;

    update_patch = patch_ttv(gather(new_img_x),Data.llr,para)*0.5;
    update_bins = compute_tTV_bins_no_lr(gather(new_img_x),weight_tTV,beta_sqrd,para.Recon.bins)*0.5;
    update_bins(permute(Data.llr.mask,[1,2,4,3])) = update_patch(permute(Data.llr.mask,[1,2,4,3]));
    update_term = update_term + update_bins; clear update_bins update_patch
    
    update_term = update_term + (Data.first_guess - new_img_x - Data.Y)*lr_weight;

%%%%% conjugate gradient
    tic;
    if iter_no > 1
        beta = update_term(:)'*update_term(:)/(update_term_old(:)'*update_term_old(:)+eps('single'));
        update_term = update_term + beta*update_term_old;
    end
    update_term_old = update_term; clear update_term
    
%%%%% line search    
    
    para.Cost = Cost_STCR(fidelity_norm, new_img_x, weight_sTV, weight_tTV, para.Cost); clear fidelity_update
    step_size = line_search_iso(new_img_x,update_term_old,Data,para);
    para.Recon.step_size(iter_no) = step_size;

    new_img_x = new_img_x + step_size * update_term_old;
    para.CPUtime.update(iter_no) = toc; 

%%%%% plot&save part 

    if ifplot == 1
        showImage3D(new_img_x,para.Cost)
    end

%%%%% break when step size too small or cost not changing too much

    if para.Recon.break && iter_no > 1
        if step_size<1e-4 %|| abs(para.Cost.totalCost(end) - para.Cost.totalCost(end-1))/para.Cost.totalCost(end-1) < 1e-4
            break
        end
    end
    
%     if mod(iter_no,1) == 0
        fprintf(['Iteration = ' num2str(iter_no) '...']);
        toc(t1);
%     end
end

Image = squeeze(gather(new_img_x));
%para = get_CPU_time(para);
%fprintf(['Iterative STCR running time is ' num2str(para.CPUtime.interative_recon) 's' '\n'])
end

