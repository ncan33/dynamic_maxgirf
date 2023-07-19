function [Image,para] = iteration_recon_3D(first_est,kSpace,sens_map,para)

disp('Performing iterative STCR reconstruction...');
disp('Showing iteration number...')

save_frequency = para.setting.save_frequency; 
ifBLS          = para.setting.BacktrackingLineSearch;
ifplot         = para.setting.ifplot;
ifGPU          = para.setting.ifGPU;
debug          = para.setting.debug;
weight_tTV     = para.weight_tTV;
weight_sTV     = para.weight_sTV;
weight_slTV    = para.weight_sliceTV;
beta_sqrd      = para.Recon.epsilon;
noi_start      = para.Recon.noi_start;
noi_end        = para.Recon.noi_end;
step_size      = para.Recon.step_size;
lambda_vbm4d   = para.vbm4d;
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

new_img_x = single(first_est);
new_img_z = single(first_est);
clear first_est

t = single(zeros(noi_end+1,1));
t(1) = 1;
for i=1:noi_end
    t(i+1) = 0.5*(1+sqrt(1+4*t(i)));
end

%sens_map_conj = conj(sens_map);
mask = logical(abs(kSpace(:,:,:,:,1)));

%%%
para.Recon.kSpace_points = sum(mask(:))*8;
para.Recon.kSpace_full_points = numel(kSpace);
para.Recon.undersampling_ratio = para.Recon.kSpace_points/para.Recon.kSpace_full_points;
para.Recon.kSpace_sum = sum(abs(kSpace(:)));
k_temp = kSpace(:);
k_temp(k_temp == 0) = [];
para.Recon.kSpace_std = std(k_temp);
para.Recon.kSpace_mean = mean(abs(k_temp));
i_temp = ifft3(kSpace);
para.Recon.ImSpace_mean = mean(abs(i_temp(:)));
para.Recon.ImSpace_std = std(abs(i_temp(:)));
clear k_temp i_temp
%%%

if ifGPU
    new_img_x = gpuArray(new_img_x);
    sens_map = gpuArray(sens_map);
    kSpace = gpuArray(kSpace);
    %sens_map_conj = gpuArray(sens_map_conj);
    mask = gpuArray(mask);
    new_img_z = gpuArray(new_img_z);
    beta_sqrd = gpuArray(beta_sqrd);
end
%%%
%ImkSpace = sum(bsxfun(@times,ifft3(kSpace),conj(sens_map)),para.Recon.kSpace_dim);
para.Cost.fidelityNorm = [];
para.Cost.temporalNorm = [];
para.Cost.spatialNorm = [];
para.Cost.totalCost = [];
%%%
for iter_no = noi_start+1:noi_end

    t1=tic;
    fprintf('%.2f%%...',iter_no/noi_end*100);

%%%%% fidelity term/temporal/spatial TV

    tic; 
    kSpace_update = bsxfun(@times,new_img_x,sens_map);
    kSpace_update = fft3(kSpace_update);
    fidelity_update = kSpace - kSpace_update; clear kSpace_update;
    fidelity_update = bsxfun(@times,fidelity_update,mask);
    fidelity_update = ifft3(fidelity_update);
    fidelity_update = sum(bsxfun(@times,fidelity_update,conj(sens_map)),para.Recon.kSpace_dim);
    para.CPUtime.fidelity(iter_no) = toc;

%{
    tic;
    fidelity_update = bsxfun(@times,new_img_x,sens_map);
    fidelity_update = fft3(fidelity_update);
    fidelity_update = bsxfun(@times,fidelity_update,mask);
    fidelity_update = ifft3(fidelity_update);
    fidelity_update = bsxfun(@times,fidelity_update,conj(sens_map));
    fidelity_update = sum(fidelity_update,para.Recon.kSpace_dim);
    fidelity_update = ImkSpace - fidelity_update; %this is like doing a denoise! get rid of k-space and using only the image fomr RAW kspace!
    %fidelity_update = ImkSpace - new_img_x;
    para.CPUtime.fidelity(iter_no) = toc;
  %}   
    tic; tTV_update = compute_3DtTV_yt(new_img_x,weight_tTV,beta_sqrd);  para.CPUtime.tTV(iter_no) = toc;
    tic; sliceTV_update = compute_sliceTV_yt(new_img_x,weight_slTV,beta_sqrd);  %para.CPUtime.sTV(iter_no) = toc;
    sTV_update = compute_sTV_yt(new_img_x,weight_sTV,beta_sqrd);  para.CPUtime.sTV(iter_no) = toc;

    %tic; sTV_update = compute_sTV_3D_yt(new_img_x,weight_sTV,beta_sqrd); para.CPUtime.sTV(iter_no) = toc;
     
    tic;
    if debug        
        para.Cost = NormCalculation3D(fidelity_update, new_img_x, weight_sTV, weight_tTV, para.Cost);
        if iter_no>1 && ifBLS==1
            if para.Cost.totalCost(iter_no)>para.Cost.totalCost(iter_no-1)
                step_size = step_size*0.5;
                %disp('Chaning step size into');disp(step_size);
            end
        end 
    end
    update_term  = fidelity_update + tTV_update + sTV_update + sliceTV_update; clear fidelity_update tTV_update sTV_update sliceTV_update

%%%%% Fast gradient descent part
    old_img_z    = new_img_z;
    new_img_z    = new_img_x + step_size * update_term; clear update_term
    new_img_x    = (1+(t(iter_no)-1)/t(iter_no+1))*new_img_z + (1-t(iter_no))/t(iter_no+1)*old_img_z;
    para.CPUtime.update(iter_no) = toc;

%%%%% VBM4D
    
    tic;
    if(lambda_vbm4d~=0)
        vbm4d_term_update = VBM4D_ye(new_img_x);
        new_img_x = new_img_x + lambda_vbm4d*(vbm4d_term_update-new_img_x);
    end
    para.CPUtime.VBM4D(iter_no) = toc;
    para.Recon.step_size(iter_no) = step_size;

%%%%% plot&save part 

    if ifplot
        showImage3D(new_img_x,para.Cost)
    end
    
    if mod(iter_no,save_frequency)==0
        Image = squeeze(new_img_x);
        disp('Saving image into Results...')
        save_dir_temp = strcat(para.Recon.save_dir,num2str(iter_no,'N%04.f_'),para.time);
        save([save_dir_temp,'.mat'],'Image','new_img_z','para','-v7.3');
    end
    toc(t1);
end

Image = squeeze(new_img_x);
para = get_CPU_time(para);
fprintf(['Iterative STCR running time is ' num2str(para.CPUtime.interative_recon) 's' '\n'])
disp('Reconstruction done');fprintf('\n')
disp('Saving image into Results...')
Image = abs(Image);
Image = crop_half_FOV(Image);
Image = gather(Image);

save_dir = strcat(para.dir.save_recon_img_mat_dir,para.dir.save_recon_img_name);
save([save_dir,'.mat'],'Image','para','-v7.3');

disp('All done.');fprintf('\n')