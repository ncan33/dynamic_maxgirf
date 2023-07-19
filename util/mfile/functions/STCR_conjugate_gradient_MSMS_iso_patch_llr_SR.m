function [Image,para] = STCR_conjugate_gradient_MSMS_iso_patch_llr_SR(Data,Data_off,para)
%[Image,para] = STCR_conjugate_gradient(Data,para)
disp('Performing iterative STCR reconstruction...');
disp('Showing progress...')

ifplot         = para.setting.ifplot;
ifGPU          = para.setting.ifGPU;
weight_tTV     = para.Recon.weight_tTV;
weight_sTV     = para.Recon.weight_sTV;
beta_sqrd      = para.Recon.epsilon;
para.Recon.step_size = para.Recon.step_size(1);

%{
%lambda_vbm4d   = para.vbm4d;
%lambda_new    = para.recon.low_rank;
%bloc_x = 32;
%bloc_y = 32;
%tau = (0.01/5)*bloc_x;
%lambda_new = para.Recon.low_rank;
% VBM4D params

% intensity range of the video
%profile = 'np';      % V-BM4D parameter profile
%  'lc' --> low complexity
%  'np' --> normal profile
%do_wiener = 1;       % Wiener filtering
%   1 --> enable Wiener filtering
%   0 --> disable Wiener filtering
%sharpen = 1;         % Sharpening
%   1 --> disable sharpening
%  >1 --> enable sharpening
%deflicker = 1;       % Deflickering
%   1 --> disable deflickering
%  <1 --> enable deflickering
%verbose = 0;         % Verbose mode
%p_sigma     = 1;      % Noise standard deviation. it should be in the same
%}
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

para.Recon.nset = length(Data);
order = vec([1:para.Recon.nSMS:para.Recon.nSMS*para.Recon.nset]' + [0,para.Recon.nset-1:-1:1]);
[~,order_back] = sort(order); 
order_back = reshape(order_back,[para.Recon.nSMS,para.Recon.nset]);


siz = size(Data{1}.first_est);
siz(5) = para.Recon.nset*para.Recon.nSMS*2+1;
new_img_x = zeros(siz,'like',Data{1}.first_est);
nslice_low_res = para.Recon.nset*para.Recon.nSMS;


para.Recon.order_back = order_back;
para.Recon.nslice_low_res = nslice_low_res;
para.Recon.siz = siz;

% if isfield(Data,'first_guess')
%     new_img_x = Data.first_guess;   
% else
%     new_img_x = single(Data.first_est);
% end

if isfield(Data{1},'phase_mod')
    Data.phase_mod_conj = conj(single(Data.phase_mod));
end

if isfield(Data{1},'sens_map')
    for i=1:para.Recon.nset
        Data{i}.sens_map_conj = conj(Data{i}.sens_map);
        Data_off{i}.sens_map_conj = conj(Data_off{i}.sens_map);
    end
end

if ifGPU
%    Data.kSpace        = gpuArray(Data.kSpace);
    new_img_x          = gpuArray(new_img_x);
    for i=1:para.Recon.nset
        Data{i}.sens_map      = gpuArray(Data{i}.sens_map);
        Data{i}.sens_map_conj = gpuArray(Data{i}.sens_map_conj);
        Data_off{i}.sens_map      = gpuArray(Data_off{i}.sens_map);
        Data_off{i}.sens_map_conj = gpuArray(Data_off{i}.sens_map_conj);
    end
%     if isfield(Data,'mask')
%         Data.mask          = gpuArray(Data.mask);
%     end
    if isfield(Data{1},'filter')
        for i=1:para.Recon.nset
            Data{i}.filter        = gpuArray(Data{i}.filter);
            Data_off{i}.filter        = gpuArray(Data_off{i}.filter);
        end
    end
    %Data.first_est = gpuArray(Data.first_est);
%    Data.phase_mod = gpuArray(Data.phase_mod);
%    Data.phase_mod_conj = gpuArray(Data.phase_mod_conj);
%    beta_sqrd = gpuArray(beta_sqrd);
    if isfield(Data{1},'N')
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

% fidelity = @(im) compute_fidelity_yt_new(im,Data,para);
spatial  = @(im) compute_sTV_3D_iso(im,weight_sTV,beta_sqrd);
temporal = @(im) compute_tTV_yt(im,weight_tTV,beta_sqrd);

for iter_no = 1:para.Recon.noi

    if mod(iter_no,10) == 1
        t1 = tic;
    end

%%%%% fidelity term/temporal/spatial TV
%%%%% super resolution fidelity
    Image_1 = reshape(new_img_x(:,:,:,:,1:end-1),[siz(1:4),2,nslice_low_res]);
    Image_1 = permute(sum(Image_1,5),[1,2,3,4,6,5]);
    Image_2 = reshape(new_img_x(:,:,:,:,2:end),[siz(1:4),2,nslice_low_res]);
    Image_2 = permute(sum(Image_2,5),[1,2,3,4,6,5]);
    
    fidelity_update_1 = zeros(size(Image_1),'like',Image_1);
    fidelity_update_2 = zeros(size(Image_2),'like',Image_2);
    for i=1:para.Recon.nset
        [fidelity_update_1(:,:,:,:,order_back(:,i)),fidelity_norm_1(i)] = compute_fidelity_yt_new(Image_1(:,:,:,:,order_back(:,i)),Data{i},para);
        [fidelity_update_2(:,:,:,:,order_back(:,i)),fidelity_norm_2(i)] = compute_fidelity_yt_new(Image_2(:,:,:,:,order_back(:,i)),Data_off{i},para);
    end
    
    fidelity_update_1 = repmat(fidelity_update_1,[1,1,1,1,1,2]);
    fidelity_update_1 = permute(fidelity_update_1,[1,2,3,4,6,5]);
    fidelity_update_1 = reshape(fidelity_update_1,[siz(1:4),siz(5)-1]);
    fidelity_update_1 = cat(5,fidelity_update_1,zeros(siz(1:4)))/4;
    
    fidelity_update_2 = repmat(fidelity_update_2,[1,1,1,1,1,2]);
    fidelity_update_2 = permute(fidelity_update_2,[1,2,3,4,6,5]);
    fidelity_update_2 = reshape(fidelity_update_2,[siz(1:4),siz(5)-1]);
    fidelity_update_2 = cat(5,zeros(siz(1:4)),fidelity_update_2)/4;
    
    update_term = fidelity_update_1 + fidelity_update_2;
    fidelity_norm = sos([fidelity_norm_1,fidelity_norm_2])/10^2;
    clear fidelity_update_1 fidelity_update_2
%%%% 
    
%     [update_term,fidelity_norm] = fidelity(iso_forward(new_img_x,Data.iso));
%     update_term = iso_backward(update_term,Data.iso);
    update_term = update_term + temporal(new_img_x);
    update_term = update_term + permute(spatial(squeeze(new_img_x)),[1,2,3,5,4]);

    %update_term = update_term + patch_low_rank(new_img_x,Data)*0.25;
%     if isfield(para.Recon,'bins')
%         update_patch = patch_low_rank(new_img_x,Data,para)*0.25;
%         update_bins = compute_tTV_bins(new_img_x,weight_tTV,beta_sqrd,para.Recon.bins)*0.25;
%         update_bins(permute(Data.llr.mask,[1,2,4,3])) = update_patch(permute(Data.llr.mask,[1,2,4,3]));
%         update_term = update_term + update_bins*2;
%         %update_term = update_term + compute_tTV_bins(new_img_x,weight_tTV,beta_sqrd,para.Recon.bins)*0.25;
%     end

%%%%% conjugate gradient
    tic;
    if iter_no > 1
        beta = update_term(:)'*update_term(:)/(update_term_old(:)'*update_term_old(:)+eps('single'));
        update_term = update_term + beta*update_term_old;
    end
    update_term_old = update_term; clear update_term
    
%%%%% line search    
    
    para.Cost = Cost_STCR(fidelity_norm, new_img_x, weight_sTV, weight_tTV, para.Cost); clear fidelity_update
    step_size = line_search_sr(new_img_x,update_term_old,Data,Data_off,para);
    para.Recon.step_size(iter_no) = step_size;

    new_img_x = new_img_x + step_size * update_term_old;
    para.CPUtime.update(iter_no) = toc;
%{
%%%%% VBM4D

    tic;
    if(lambda_vbm4d~=0)
        vbm4d_term_update = VBM4D_ye(new_img_x);
        new_img_x = new_img_x + lambda_vbm4d*(vbm4d_term_update-new_img_x);
    end
    para.CPUtime.VBM4D(iter_no) = toc;
 
%%%%% Low Rank

    if lambda_new~=0
        for i=1:size(new_img_x,5)
            lr_update(:,:,:,:,i) = low_rank_yt(new_img_x(:,:,:,:,i),bloc_x,bloc_y,tau);
        end

        new_img_x = new_img_x + lambda_new * (lr_update - new_img_x);
        clear lr_update
    end
%}
%%%%% plot&save part 

    if ifplot == 1
        showImage(new_img_x,para.Cost)
    end

%%%%% break when step size too small or cost not changing too much

    if para.Recon.break && iter_no > 1
        if step_size<1e-4 %|| abs(para.Cost.totalCost(end) - para.Cost.totalCost(end-1))/para.Cost.totalCost(end-1) < 1e-4
            break
        end
    end
    
    if mod(iter_no,10) == 0
        fprintf(['Iteration = ' num2str(iter_no) '...']);
        toc(t1);
    end
end

Image = squeeze(gather(new_img_x));
%para = get_CPU_time(para);
%fprintf(['Iterative STCR running time is ' num2str(para.CPUtime.interative_recon) 's' '\n'])