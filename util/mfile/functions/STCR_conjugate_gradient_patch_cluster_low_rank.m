function [Image,para] = STCR_conjugate_gradient_patch_cluster_low_rank(Data,para)
% [Image,para] = iteration_recon(first_est,kSpace,sens_map,phase_mod,para,varargin)
disp('Performing iterative STCR reconstruction...');
disp('Showing progress...')

% save_frequency = para.setting.save_frequency;
%ifNUFFT        = para.ifNUFFT; 
%ifBLS          = para.setting.BacktrackingLineSearch;
ifplot         = para.setting.ifplot;
ifGPU          = para.setting.ifGPU;
%debug          = para.setting.debug;
%nSMS           = para.Recon.nSMS;
weight_tTV     = para.Recon.weight_tTV;
weight_sTV     = para.Recon.weight_sTV;
%weight_slTV    = para.Recon.weight_sTV/2;
%para.Recon.weight_slTV = weight_slTV;
beta_sqrd      = para.Recon.epsilon;
% noi_start      = para.Recon.noi_start;
% noi_end        = para.Recon.noi_end;
para.Recon.step_size = para.Recon.step_size(1);
% lambda_vbm4d   = para.vbm4d;
%lambda_new    = para.recon.low_rank;

%
bloc_x = 8;
bloc_y = 8;
% tau = (0.01/5)*bloc_x;
tau = max(abs(Data.first_est(:)))*0.5;
tau = 2000;
lambda_new = 0.125;

if isfield(para.Recon,'RF_frames')
    RF_frames = para.Recon.RF_frames;
    PD_frames = para.Recon.PD_frames;
    Data.first_est = Data.first_est(:,:,RF_frames,:,:);
    if isfield(Data,'mask')
        Data.mask = Data.mask(:,:,RF_frames,:,:,:,:);
    end
    if isfield(Data,'phase_mod')
        Data.phase_mod = Data.phase_mod(:,:,RF_frames,:,:,:,:);
    end
    if isfield(Data,'N')
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

%Data.first_est = Data.first_est/max(abs(Data.first_est(:)));
if not(para.Recon.ifContinue)
    new_img_x = single(Data.first_est);
    %new_img_x = repmat(mean(single(Data.first_est),3),[1 1 67]);
    
    %new_img_x = repmat(single(Data.first_est),[1 1 1 1 1 5]);
    %new_img_x(:,:,:,1) = circshift(new_img_x(:,:,:,1),-2,3);
    %new_img_x(:,:,:,2) = circshift(new_img_x(:,:,:,1),-1,3);
    %new_img_x(:,:,:,3) = circshift(new_img_x(:,:,:,1),1,3);
    %new_img_x(:,:,:,4) = circshift(new_img_x(:,:,:,1),2,3);
    %new_img_x = mean(new_img_x,6);
    
else
    load(para.dir.continue_dir,'Image','fidelityNorm','spatialNorm','temporalNorm','totalCost')
    new_img_x = permute(Image,[1 2 3 5 4]); clear Image
end

if isfield(Data,'first_guess')
    new_img_x = Data.first_guess;
end

%t = single(zeros(noi_end+1,1));
%t(1) = 1;
%for i=1:noi_end
%    t(i+1) = 0.5*(1+sqrt(1+4*t(i)));
%end

if isfield(Data,'phase_mod')
    Data.phase_mod_conj = conj(single(Data.phase_mod));
end
if isfield(Data,'sens_map')
    Data.sens_map_conj = conj(Data.sens_map);
end

%{
k_temp = kSpace(:);
para.Recon.kSpace_points = sum(logical(abs(k_temp)));
para.Recon.kSpace_full_points = numel(kSpace);
para.Recon.kSpace_sum = sum(abs(kSpace(:)));
para.Recon.undersampling_ratio = para.Recon.kSpace_points/para.Recon.kSpace_full_points;
k_temp(k_temp == 0) = [];
para.Recon.kSpace_std = std(k_temp);
para.Recon.kSpace_mean = mean(abs(k_temp));
i_temp = ifft2(kSpace);
para.Recon.ImSpace_mean = mean(abs(i_temp(:)));
para.Recon.ImSpace_std = std(abs(i_temp(:)));
clear k_temp i_temp
%}
%gpuDevice(3);
if ifGPU
    new_img_x = gpuArray(new_img_x);
    Data.first_est = gpuArray(Data.first_est);
    %Data.sens_map = gpuArray(Data.sens_map);
    if isfield(Data,'N')
        for i=1:length(Data.N)
            Data.N(i).S = gpuArray(Data.N(i).S);
            Data.N(i).Apodizer = gpuArray(Data.N(i).Apodizer);
            Data.N(i).W = gpuArray(Data.N(i).W);
        end
    end
    %img_k = gpuArray(img_k);
    %kSpace = gpuArray(kSpace);
    %Data.phase_mod_conj = gpuArray(Data.phase_mod_conj);
    %Data.sens_map_conj = gpuArray(Data.sens_map_conj);
    %new_img_z = gpuArray(new_img_z);
    %beta_sqrd = gpuArray(beta_sqrd);
end
temporal = @(im) compute_tTV_yt(im,weight_tTV,beta_sqrd);

para.Cost = struct('fidelityNorm',[],'temporalNorm',[],'spatialNorm',[],'totalCost',[]);

for iter_no = 1:para.Recon.noi
   
    t1 = tic;
    %fprintf('%.2f%%...',iter_no/noi_end*100);
    fprintf(['Iteration = ' num2str(iter_no) '...']);

%%%%% fidelity term/temporal/spatial TV

    tic; 
    [update_term,fidelity_norm] = compute_fidelity_yt_new(new_img_x,Data,para);
    para.CPUtime.fidelity(iter_no) = toc;
    
    tic; 
    update_term = update_term + temporal(new_img_x) * 0.5;
    update_term = update_term + weight_tTV * apply_patch_cluster_ttv(new_img_x, Data.patch) * 0.5;
    para.CPUtime.tTV(iter_no) = toc;
    
    tic; 
    update_term = update_term + compute_sTV_yt(new_img_x,weight_sTV,beta_sqrd);
    para.CPUtime.sTV(iter_no) = toc;

    update_term = update_term + lambda_new * (apply_patch_cluster_low_rank(new_img_x, Data.patch, tau) - new_img_x);


%%%%% conjugate gradient
    tic;
    if iter_no > 1
        beta = update_term(:)'*update_term(:)/(update_term_old(:)'*update_term_old(:)+eps('single'));
        update_term = update_term + beta*update_term_old;
    end
    update_term_old = update_term; clear update_term
    
%%%%% line search    
    
    para.Cost = Cost_STCR_cluster(fidelity_norm, new_img_x, weight_sTV, weight_tTV, Data.patch, para.Cost); clear fidelity_update
    step_size = line_search(new_img_x,update_term_old,Data,para);
    para.Recon.step_size(iter_no) = step_size;

    new_img_x   = new_img_x + step_size * update_term_old;
    para.CPUtime.update(iter_no) = toc;

    
%     new_img_x = new_img_x + lambda_new * (update_lr - new_img_x);
    clear update_lr

%%%%% VBM4D
    
%     tic;
%     if(lambda_vbm4d~=0)
%         vbm4d_term_update = VBM4D_ye(new_img_x);
%         new_img_x = new_img_x + lambda_vbm4d*(vbm4d_term_update-new_img_x);
%     end
%     para.CPUtime.VBM4D(iter_no) = toc;

%%%%% plot&save part 

    if ifplot ==1
        showImage(new_img_x,para.Cost)
    end
    
    % break when step size too small or cost not changing too much
    if iter_no > 1
        if step_size<1e-4 %|| abs(para.Cost.totalCost(end) - para.Cost.totalCost(end-1))/para.Cost.totalCost(end-1) < 1e-4
            break
        end
    end
    
%     if mod(iter_no,save_frequency)==0
%         Image = squeeze(new_img_x);
%         disp('Saving image into Results...')
%         save_dir_temp = strcat(para.Recon.save_dir,para.time);
%         save([save_dir_temp,'.mat'],'Image','new_img_z','para','-v7.3');
%     end
    
    toc(t1);
end

Image = squeeze(gather(new_img_x));
para = get_CPU_time(para);
fprintf(['Iterative STCR running time is ' num2str(para.CPUtime.interative_recon) 's' '\n'])
end



function [Cost_new,Cost,fNorm,tNorm,sNorm] = Cost_STCR_cluster(fUpdate, Image, sWeight, tWeight, patch, Cost_old)

N = numel(Image);

fNorm = sum(abs(fUpdate(:)).^2);
% Image = crop_half_FOV(Image);
if tWeight ~= 0
    tNorm = 0.5 * mean(tWeight(:)) .* abs(diff(Image,1,3));
    tNorm = sum(tNorm(:));
else
    tNorm = 0;
end

patch_size = size(patch.idx{1});
n_patch = size(patch.idx);
n_cluster = size(patch.cluster, 2);
ttv = zeros(size(Image), 'like', Image);
for ipatch = 1:n_patch
    patch_temp = Image(patch.idx{ipatch});
    for icluster = 1:n_cluster
        if length(patch.cluster{ipatch, icluster}) > 3
            ttv(patch.idx{ipatch}(:, :, patch.cluster{ipatch, icluster})) = ttv(patch.idx{ipatch}(:, :, patch.cluster{ipatch, icluster})) + cat(3, abs(diff(patch_temp(:,:,patch.cluster{ipatch, icluster}), 1, 3)), zeros(patch_size(1:2), 'like', Image));
        end
    end
end
ttv = ttv .* patch.mask;
tNorm = tNorm + sum(ttv(:)) * tWeight * 0.5;

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

fNorm = fNorm/N;
tNorm = tNorm/N;
sNorm = sNorm/N;

Cost = sNorm + tNorm + fNorm;

if nargin == 5
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


function step = line_search(old, update, Data, para)
%--------------------------------------------------------------------------
%   [step] = line_search(old, update, Data, para)
%--------------------------------------------------------------------------
%   Line search called in a conjugate gradient algorithm
%--------------------------------------------------------------------------
%   Inputs:      
%       - old       [sx, sy, nof, ...]
%       - update    [sx, sy, nof, ...]
%       - Data      [structure]
%       - para      [structure]
%               
%       - old       image from previous iteration
%       - update    update term
%       - Data      see 'help STCR_conjugate_gradient.m'
%       - para      see 'help STCR_conjugate_gradient.m'
%--------------------------------------------------------------------------
%   Output:
%       - step      [scalar]
%
%       - step      step size for CG update
%--------------------------------------------------------------------------
%   This function trys to find a suitable step size to perform a CG update.
%   The function starts with a step size adopted from last iteration, and
%   multiply it by 1.3 (magic number). If the step size yeilds a cost that
%   is larger than the previous cost, it shrinks the step size by 0.8
%   (magic number again). If it yeilds a cost that is smaller than the
%   previous cost, it will increase the step size by 1.3 until it no longer
%   yeild a smaller cost. The maximum number of trys is 15.
%--------------------------------------------------------------------------
%   Author:
%       Ye Tian
%       E-mail: phye1988@gmail.com
%--------------------------------------------------------------------------

step_start = para.Recon.step_size(end)*1.3; % magic number
%step_start = 2;
%step_start = para.Recon.step_size(1);
tau = 0.8; % magic number again
tau_2 = 1.3;
max_try = 20;
step = step_start;

cost_old = para.Cost.totalCost(end);
flag = 0;

for i=1:max_try
%      fprintf(['Iter = ' num2str(i) '... '])
    
    new = old + step * update;
    fidelity_new = compute_fidelity_for_line_search_yt(new, Data, para);
    cost_new = Cost_STCR_cluster(fidelity_new, new, para.Recon.weight_sTV, para.Recon.weight_tTV, Data.patch);

%     fprintf(['Cost new = ' num2str(round(cost_new)) '...\n'])
    if cost_new > cost_old && flag == 0
        step = step * tau;
    elseif cost_new < cost_old 
        step = step * tau_2;
        cost_old = cost_new;
        flag = 1;
    elseif cost_new > cost_old && flag == 1
        step = step / tau_2;
%          fprintf(['Step = ' num2str(step) '...\n'])
%          fprintf(['Cost = ' num2str(round(cost_old)) '...\n'])
        return
    end
end
%  fprintf(['Step = ' num2str(step) '...\n'])
%  fprintf(['Cost = ' num2str(round(cost_new)) '...\n'])
end