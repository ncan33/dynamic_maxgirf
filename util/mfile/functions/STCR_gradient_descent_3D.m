function [Image,para] = STCR_gradient_descent_3D(Data,para)
% [Image,para] = iteration_recon(first_est,kSpace,sens_map,phase_mod,para,varargin)
disp('Performing iterative STCR reconstruction...');
disp('Showing progress...')

% save_frequency = para.setting.save_frequency;
%ifNUFFT        = para.ifNUFFT; 
% ifBLS          = para.setting.BacktrackingLineSearch;
ifplot         = para.setting.ifplot;
ifGPU          = para.setting.ifGPU;
% debug          = para.setting.debug;
%nSMS           = para.Recon.nSMS;
weight_tTV     = para.Recon.weight_tTV;
weight_sTV     = para.Recon.weight_sTV;
weight_slTV    = para.Recon.weight_sliceTV;
beta_sqrd      = para.Recon.epsilon;
%noi_start      = para.Recon.noi_start;
%noi_end        = para.Recon.noi_end;
para.Recon.step_size = para.Recon.step_size(1);
step_size      = para.Recon.step_size;
% lambda_vbm4d   = para.vbm4d;
%lambda_new    = para.recon.low_rank;

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
    %Data.first_est = gpuArray(Data.first_est);
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

para.Cost = struct('fidelityNorm',[],'temporalNorm',[],'spatialNorm',[],'totalCost',[]);
fprintf(' Iteration       Cost       Step    Time(s) \n')
for iter_no = 1:para.Recon.noi
   
    t1 = tic;
    %fprintf('%.2f%%...',iter_no/noi_end*100);
%     fprintf(['Iteration = ' num2str(iter_no) '...']);

%%%%% fidelity term/temporal/spatial TV

    tic; fidelity_update = compute_fidelity_yt_new(new_img_x,Data,para); para.CPUtime.fidelity(iter_no) = toc;
    tic; tTV_update = compute_3DtTV_yt(new_img_x,weight_tTV,beta_sqrd); para.CPUtime.tTV(iter_no) = toc;
    tic; sTV_update = compute_sTV_yt(new_img_x,weight_sTV,beta_sqrd);
    slTV_update = compute_sliceTV_yt(new_img_x,weight_slTV,beta_sqrd);  para.CPUtime.sTV(iter_no) = toc;
    
    update_term = fidelity_update + tTV_update + sTV_update + slTV_update;
    clear tTV_update sTV_update slTV_update


    tic;    
%     if debug
         para.Cost = Cost_STCR_3D(fidelity_update, new_img_x, weight_sTV, weight_tTV, weight_slTV, para.Cost); clear fidelity_update
        %if iter_no>1 && abs(para.Cost.totalCost(end) - para.Cost.totalCost(end-1))/para.Cost.totalCost(end-1) < 1e-3
            %para.Recon.step_size(iter_no) = 0.02;
        %end
%         fprintf([num2str(round(para.Cost.totalCost(end))) '...'])

%        if iter_no>1 && ifBLS==1
%            if para.Cost.totalCost(iter_no)>para.Cost.totalCost(iter_no-1)
%                step_size = step_size*0.8;
%            else
%                step_size = step_size*1.1;
%            end
%        end 
%        step_size = line_search_golden_section(new_img_x,update_term,Data,para);
%         para.Recon.step_size(iter_no) = step_size;
%         fprintf(num2str(step_size))
%     end

    new_img_x   = new_img_x + step_size * update_term;
    para.CPUtime.update(iter_no) = toc;

%%%%% VBM4D
    
%    tic;
%    if(lambda_vbm4d~=0)
%        vbm4d_term_update = VBM4D_ye(new_img_x);
%        new_img_x = new_img_x + lambda_vbm4d*(vbm4d_term_update-new_img_x);
%    end
%    para.CPUtime.VBM4D(iter_no) = toc;

%%%%% plot&save part 

    if ifplot ==1
        showImage3D(new_img_x, para.Cost)
    end
    
    % break when step size too small or cost not changing too much
%    if iter_no > 1
%        if step_size<1e-4 || abs(para.Cost.totalCost(end) - para.Cost.totalCost(end-1))/para.Cost.totalCost(end-1) < 1e-4
%            break
%        end
%    end
    
%     if mod(iter_no,save_frequency)==0
%         Image = squeeze(new_img_x);
%         disp('Saving image into Results...')
%         save_dir_temp = strcat(para.Recon.save_dir,para.time);
%         save([save_dir_temp,'.mat'],'Image','new_img_z','para','-v7.3');
%     end
    

fprintf(sprintf('%10.0f %10.2f %10.4f %10.2f \n',iter_no,para.Cost.totalCost(end),step_size,toc(t1)));
end

Image = squeeze(new_img_x);
para = get_CPU_time(para);
fprintf(['Iterative STCR running time is ' num2str(para.CPUtime.interative_recon) 's' '\n'])
