function [Image,para] = STCR_conjugate_gradient_MSMS_simple_motion_lr(Data,para)
%[Image,para] = STCR_conjugate_gradient(Data,para)
disp('Performing iterative STCR reconstruction...');
disp('Showing progress...')

ifplot         = para.setting.ifplot;
ifGPU          = para.setting.ifGPU;
weight_tTV     = para.Recon.weight_tTV;
weight_sTV     = para.Recon.weight_sTV;
beta_sqrd      = para.Recon.epsilon;
para.Recon.step_size = para.Recon.step_size(1);




if isfield(para.Recon,'RF_frames')
    RF_frames = para.Recon.RF_frames;
    PD_frames = para.Recon.PD_frames;
    Data.first_est = Data.first_est(:,:,RF_frames,:,:);
    Data.kSpace = Data.kSpace(:,:,RF_frames,:,:,:,:);
    if isfield(Data,'mask')
        Data.mask = Data.mask(:,:,RF_frames,:,:,:,:);
    end
    if isfield(Data,'phase_mod')
        Data.phase_mod = Data.phase_mod(:,:,RF_frames,:,:,:,:);
    end
    if isfield(Data,'N')
        if PD_frames(end) == 0
            PD_frames = 1:sum(PD_frames);
        end
        Data.N.S = Data.N.S(Data.N.sx_over.^2*PD_frames(end)+1:end,Data.N.siz(1)*Data.N.siz(2)*PD_frames(end)+1:end);
        Data.N.siz(3) = size(Data.first_est,3);
%        Data.N.W = Data.N.W(:,:,RF_frames);
    end
elseif isfield(para.Recon,'PD_frames') && sum(para.Recon.PD_frames) == length(para.Recon.PD_frames)
    Image = [];
    return
end

if isfield(Data,'first_guess')
    new_img_x = Data.first_guess;   
else
    new_img_x = single(Data.first_est);
end


[sx,sy,nof,~,nSMS,nset] = size(new_img_x);
nslice = nSMS*nset;
new_img_x = reshape(new_img_x,[sx,sy,nof,nslice]);
order = 1:nSMS:nslice;
for i=nSMS:-1:2
    order = [order,i:nSMS:nslice];
end
[~,order_back] = sort(order);
iso.order = order;
iso.order_back = order_back;
iso.nSMS = nSMS;
iso.nset = nset;
Data.iso = iso;
new_img_x = new_img_x(:,:,:,order);
new_img_x = permute(new_img_x,[1,2,4,3]);
new_img_x = linear_shifts_kspace(new_img_x,Data.shifts);




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
spatial  = @(im) compute_sTV_3D_iso(im,weight_sTV,beta_sqrd);
temporal = @(im) compute_3DtTV_yt(im,weight_tTV,beta_sqrd);


for iter_no = 1:para.Recon.noi

%     if mod(iter_no,1) == 1
        t1 = tic;
%     end

%%%%% fidelity term/temporal/spatial TV

    [update_term,fidelity_norm] = fidelity(iso_forward(permute(linear_shifts_kspace(new_img_x,-Data.shifts),[1,2,4,3]),iso));
    update_term = iso_backward(update_term,iso);
    update_term = permute(update_term,[1,2,4,3]);
    update_term = linear_shifts_kspace(update_term,Data.shifts);
    update_term = update_term + temporal(new_img_x)*0.5;
    update_term = update_term + spatial(new_img_x);

    %update_term = update_term + patch_low_rank(new_img_x,Data)*0.25;
    if isfield(para.Recon,'bins')
        %update_patch = patch_low_rank(new_img_x,Data,para)*0.5;
        update_bins = compute_tTV_bins(permute(new_img_x,[1,2,4,3]),weight_tTV,beta_sqrd,para.Recon.bins)*0.5;
        %update_bins(permute(Data.llr.mask,[1,2,4,3])) = update_patch(permute(Data.llr.mask,[1,2,4,3]));
        update_term = update_term + permute(update_bins,[1,2,4,3]);
        %update_term = update_term + compute_tTV_bins(new_img_x,weight_tTV,beta_sqrd,para.Recon.bins)*0.25;
    end

%%%%% conjugate gradient
    tic;
    if iter_no > 1
        beta = update_term(:)'*update_term(:)/(update_term_old(:)'*update_term_old(:)+eps('single'));
        update_term = update_term + beta*update_term_old;
    end
    update_term_old = update_term; clear update_term
    
%%%%% line search    

    para.Cost = Cost_STCR(fidelity_norm, permute(new_img_x,[1,2,4,3]), weight_sTV, weight_tTV, para.Cost); clear fidelity_update
    step_size = line_search_iso(permute(new_img_x,[1,2,4,3]),permute(update_term_old,[1,2,4,3]),Data,para);
    para.Recon.step_size(iter_no) = step_size;

    new_img_x = new_img_x + step_size * update_term_old;
    para.CPUtime.update(iter_no) = toc;



% %%%%% plot&save part 
% 
%     if ifplot == 1
%         showImage3D(new_img_x,para.Cost)
%     end
% 
% %%%%% break when step size too small or cost not changing too much
% 
%     if para.Recon.break && iter_no > 1
%         if step_size<1e-4 %|| abs(para.Cost.totalCost(end) - para.Cost.totalCost(end-1))/para.Cost.totalCost(end-1) < 1e-4
%             break
%         end
%     end
    
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
