function [Image,para] = STCR_conjugate_gradient_pixel_bins_interp(Data,para)
% [Image,para] = iteration_recon(first_est,kSpace,sens_map,phase_mod,para,varargin)
disp('Performing iterative STCR reconstruction...');
disp('Showing progress...')

ifplot         = para.setting.ifplot;
ifGPU          = para.setting.ifGPU;
weight_tTV     = para.Recon.weight_tTV;
weight_sTV     = para.Recon.weight_sTV;
beta_sqrd      = para.Recon.epsilon;
%para.Recon.step_size = para.Recon.step_size(1);
 
%lambda_vbm4d   = para.vbm4d;
%lambda_new    = para.recon.low_rank;

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
        Data.N.W = Data.N.W(:,:,RF_frames);
    end
end

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


if isfield(Data,'first_guess')
    new_img_x = Data.first_guess;
    noi_start = length(para.Recon.step_size)+1;
    para.Recon.noi = para.Recon.noi + noi_start;
    if ~isfield(para,'Cost')
        para.Cost = struct('fidelityNorm',[],'temporalNorm',[],'spatialNorm',[],'totalCost',[]);
    else
        para.Cost = para.Cost;
    end
else
    new_img_x = Data.first_est;
    noi_start = 1;
    para.Cost = struct('fidelityNorm',[],'temporalNorm',[],'spatialNorm',[],'totalCost',[]);
end

if isfield(Data,'phase_mod')
    Data.phase_mod_conj = conj(single(Data.phase_mod));
end

if isfield(Data,'sens_map')
    Data.sens_map_conj = conj(Data.sens_map);
end

if ifGPU
    Data.kSpace        = gpuArray(Data.kSpace);
    new_img_x          = gpuArray(new_img_x);
    Data.sens_map      = gpuArray(Data.sens_map);
    Data.sens_map_conj = gpuArray(Data.sens_map_conj);
    Data.mask          = gpuArray(Data.mask);
    Data.filter        = gpuArray(Data.filter);
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
    
    gpuInfo = gpuDevice;
    gpuSize = gpuInfo.AvailableMemory;
    imSize  = numel(new_img_x)*8;
    if imSize*para.Recon.no_comp > gpuSize*0.3
        para.Recon.type = [para.Recon.type,' less memory'];
    end
end

%para.Cost = struct('fidelityNorm',[],'temporalNorm',[],'spatialNorm',[],'totalCost',[]);

fidelity = @(im) compute_fidelity_yt_new(im,Data,para);
spatial  = @(im) compute_sTV_yt(im,weight_sTV,beta_sqrd);
ptv_all  = @(im) compute_tTV_pixel_interp(im,weight_tTV*0.5,beta_sqrd,Data.Motion);
ptv_bin  = @(im) compute_tTV_pixel_bins_interp(im,weight_tTV*0.5,beta_sqrd,Data.Motion_bins,para.Recon.bins);
for iter_no = noi_start:para.Recon.noi

    t1 = tic;
    %fprintf('%.2f%%...',iter_no/noi_end*100);
    fprintf(['Iteration = ' num2str(iter_no) '...']);

%%%%% fidelity term/temporal/spatial TV

    tic;
    [update_term,fidelity_norm] = fidelity(new_img_x);
    para.CPUtime.fidelity(iter_no) = toc;
    
    tic;
    update_term = update_term + ptv_all(new_img_x);
    update_term = update_term + ptv_bin(new_img_x);
    %tTV_update = compute_tTV_pixel(new_img_x,weight_tTV,beta_sqrd,Data.Motion)*0.5;
    %tTV_update = tTV_update + compute_tTV_pixel_bins(new_img_x,weight_tTV,beta_sqrd,Data.Motion_bins,para.Recon.bins)*0.5;
    %tTV_update = tTV_update + compute_tTV_bins(new_img_x,weight_tTV,beta_sqrd,para.Recon.bins)*0.75;
    para.CPUtime.tTV(iter_no) = toc;
    
    
    tic;
    update_term = update_term + spatial(new_img_x);
    %sTV_update = compute_sTV_yt(new_img_x,weight_sTV,beta_sqrd);
    %sTV_update = sTV_update.*Data.sTV_mask;
    %tTV_update = tTV_update.*Data.sTV_mask;
    para.CPUtime.sTV(iter_no) = toc;

    %update_term = fidelity_update + tTV_update + sTV_update;
    %clear tTV_update sTV_update slTV_update

%%%%% conjugate gradient
    tic;
    if iter_no > noi_start
        beta = update_term(:)'*update_term(:)/(update_term_old(:)'*update_term_old(:)+eps('single'));
        update_term = update_term + beta*update_term_old;
    end
    update_term_old = update_term; clear update_term
    
%%%%% line search    

    %fidelity_update = compute_fidelity_for_line_search_yt(new_img_x,Data,para);
    
    para.Cost = Cost_STCR_pixel_bins(fidelity_norm, new_img_x, Data, para, para.Cost); clear fidelity_update
    step_size = line_search_pixel_bins(new_img_x,update_term_old,Data,para);
    para.Recon.step_size(iter_no) = step_size;

    new_img_x = new_img_x + step_size .* update_term_old;
    para.CPUtime.update(iter_no) = toc;

%%%%% VBM4D
%{
    tic;
    if(lambda_vbm4d~=0)
        vbm4d_term_update = VBM4D_ye(new_img_x);
        new_img_x = new_img_x + lambda_vbm4d*(vbm4d_term_update-new_img_x);
    end
    para.CPUtime.VBM4D(iter_no) = toc;
%}
%%%%% plot&save part 

    if ifplot == 1
        showImage(new_img_x,para.Cost)
    end
    
    % break when step size too small or cost not changing too much
    if iter_no > 1 && para.Recon.break
        if step_size<1e-4 %|| abs(para.Cost.totalCost(end) - para.Cost.totalCost(end-1))/para.Cost.totalCost(end-1) < 1e-4
            break
        end
    end
%{
    if mod(iter_no,save_frequency)==0
        Image = squeeze(new_img_x);
        disp('Saving image into Results...')
        save_dir_temp = strcat(para.Recon.save_dir,para.time);
        save([save_dir_temp,'.mat'],'Image','new_img_z','para','-v7.3');
    end
%}    
    toc(t1);
end

Image = squeeze(gather(new_img_x));
para = get_CPU_time(para);
fprintf(['Iterative STCR running time is ' num2str(para.CPUtime.interative_recon) 's' '\n'])
