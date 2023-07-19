function [Image,para] = STCR_conjugate_gradient_bins(Data,para)
% [Image,para] = iteration_recon(first_est,kSpace,sens_map,phase_mod,para,varargin)
disp('Performing iterative STCR reconstruction...');
disp('Showing progress...')

disp_freq      = 1;
ifplot         = para.setting.ifplot;
ifGPU          = para.setting.ifGPU;
weight_tTV     = para.Recon.weight_tTV;
weight_sTV     = para.Recon.weight_sTV;
beta_sqrd      = para.Recon.epsilon;
para.Recon.step_size = para.Recon.step_size(1);
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
else 
    new_img_x = Data.first_est;
end

if isfield(Data,'phase_mod')
    Data.phase_mod_conj = conj(single(Data.phase_mod));
end

if isfield(Data,'sens_map')
    Data.sens_map_conj = conj(Data.sens_map);
end

if ifGPU
    new_img_x = gpuArray(new_img_x);
    Data.first_est = gpuArray(Data.first_est);
    %kSpace = gpuArray(kSpace);
    
    %Data.sens_map = gpuArray(Data.sens_map);
    if isfield(Data,'N')
        for i=1:length(Data.N)
            Data.N(i).S = gpuArray(Data.N(i).S);
            Data.N(i).Apodizer = gpuArray(Data.N(i).Apodizer);
            Data.N(i).W = gpuArray(Data.N(i).W);
        end
    end
    %img_k = gpuArray(img_k);
    
    %Data.phase_mod_conj = gpuArray(Data.phase_mod_conj);
    %Data.sens_map_conj = gpuArray(Data.sens_map_conj);
    %new_img_z = gpuArray(new_img_z);
    %beta_sqrd = gpuArray(beta_sqrd);
end

% initialize Cost
para.Cost = struct('fidelityNorm',[],'temporalNorm',[],'spatialNorm',[],'totalCost',[]);

% initialize function handels 
fidelity = @(im) compute_fidelity_yt_new(im,Data,para);
spatial  = @(im) compute_sTV_yt(im,weight_sTV,beta_sqrd);
temporal = @(im) compute_tTV_yt(im,weight_tTV,beta_sqrd);
temp_bin = @(im) compute_tTV_bins_no_lr(im,weight_tTV,beta_sqrd, para.Recon.bins);

fprintf(' Iteration       Cost       Step    Time(s) \n')

for iter_no = 1:para.Recon.noi

    if mod(iter_no,disp_freq) == 1 || iter_no == 1 || disp_freq == 1
        t1 = tic;
    end

%% fidelity term/temporal/spatial TV

    tic;
    [update_term,fidelity_norm] = fidelity(new_img_x);
    para.CPUtime.fidelity(iter_no) = toc;
    tic;
    update_term = update_term + temporal(new_img_x);
    update_term = update_term + temp_bin(new_img_x);
    para.CPUtime.tTV(iter_no) = toc;
    
    tic;
    update_term = update_term + spatial(new_img_x);
    para.CPUtime.sTV(iter_no) = toc;

%%%%% conjugate gradient
    tic;
    if iter_no > 1
        beta = update_term(:)'*update_term(:)/(update_term_old(:)'*update_term_old(:)+eps('single'));
        update_term = update_term + beta*update_term_old;
    end
    update_term_old = update_term; clear update_term
    
%%%%% line search    
    
    para.Cost = Cost_STCR_bins(fidelity_norm, new_img_x, para, para.Cost); clear fidelity_update
    step_size = line_search_bins(new_img_x,update_term_old,Data,para);
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
        if step_size<1e-4 || abs(para.Cost.totalCost(end) - para.Cost.totalCost(end-1))/para.Cost.totalCost(end-1) < 1e-4
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
    if mod(iter_no,disp_freq) == 0 || iter_no == para.Recon.noi
        fprintf(sprintf('%10.0f %10.2f %10.4f %10.2f \n',iter_no,para.Cost.totalCost(end),step_size,toc(t1)));
    end
end

Image = squeeze(gather(new_img_x));
para = get_CPU_time(para);
fprintf(['Iterative STCR running time is ' num2str(para.CPUtime.interative_recon) 's' '\n'])
