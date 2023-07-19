function [Image,para] = iteration_recon_GNUFFT(first_est,G,kSpace,kSpace_radial,sens_map,para,varargin)

disp('Performing iterative STCR reconstruction...');
disp('Showing progress...')

save_frequency = para.save_frequency;
ifNUFFT        = para.ifNUFFT; 
ifBLS          = para.BacktrackingLineSearch;
ifplot         = para.plot;
ifGPU          = para.ifGPU;
debug          = para.debug;
nSMS           = para.Recon.nSMS;
weight_tTV     = para.weight_tTV;
weight_sTV     = para.weight_sTV;
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

if not(para.ifContinue)
    new_img_x = single(first_est);
    new_img_z = single(first_est);
else
    load(para.dir.continue_dir,'Image','new_img_z','fidelityNorm','spatialNorm','temporalNorm','totalCost')
    new_img_x = permute(Image,[1 2 3 5 4]); clear Image
end
clear first_est

t = single(zeros(noi_end+1,1));
t(1) = 1;
for i=1:noi_end
    t(i+1) = 0.5*(1+sqrt(1+4*t(i)));
end

%phase_mod_conj = conj(phase_mod);
sens_map_conj = conj(sens_map);

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

if ifGPU
    new_img_x = gpuArray(new_img_x);
    sens_map = gpuArray(sens_map);
    %img_k = gpuArray(img_k);
    kSpace = gpuArray(kSpace);
%    phase_mod_conj = gpuArray(phase_mod_conj);
    sens_map_conj = gpuArray(sens_map_conj);
    new_img_z = gpuArray(new_img_z);
    beta_sqrd = gpuArray(beta_sqrd);
end

mask = logical(abs(kSpace));clear kSpace

%sx = G{1}.siz(1);
%fftshift_mask = single(ones(sx+1));
%fftshift_mask(1:2:sx+1,:) = -fftshift_mask(1:2:sx+1,:);
%fftshift_mask(:,1:2:sx+1) = -fftshift_mask(:,1:2:sx+1);

for iter_no = noi_start+1:noi_end

    t1 = tic;
    fprintf('%.2f%%...',iter_no/noi_end*100);

%%%%% fidelity term/temporal/spatial TV
   
    tic;
    kSpace_update = bsxfun(@times,new_img_x,sens_map);
    kSpace_update_radial = GROG.GNUFFT_SMS(kSpace_update,mask,G);
    fidelity_update_radial = kSpace_radial - kSpace_update_radial;
    fidelity_update = GROG.GNUFFT_adj_SMS(fidelity_update_radial,G);
 %{   
    %kSpace_update = fftshift(kSpace_update,1);
    %kSpace_update = fftshift(kSpace_update,2);
    %kSpace_update = bsxfun(@times,kSpace_update,fftshift_mask);
    %kSpace_update = fft2(kSpace_update);
    %kSpace_update = bsxfun(@times,kSpace_update,fftshift_mask);
    %kSpace_update = fftshift(kSpace_update,1);
    %kSpace_update = fftshift(kSpace_update,2);
    %kSpace_update = bsxfun(@times,kSpace_update,mask);

    kSpace_update_temp(:,:,:,:,1,1,1) = sum(kSpace_update,5);
    kSpace_update_temp(:,:,:,:,1,1,2) = kSpace_update(:,:,:,:,1,:,:) + exp(-1i*2*pi/3)*kSpace_update(:,:,:,:,2,:,:) + exp(-1i*4*pi/3)*kSpace_update(:,:,:,:,3,:,:);
    kSpace_update_temp(:,:,:,:,1,1,3) = kSpace_update(:,:,:,:,1,:,:) + exp(-1i*4*pi/3)*kSpace_update(:,:,:,:,2,:,:) + exp(-1i*2*pi/3)*kSpace_update(:,:,:,:,3,:,:);
    kSpace_update_temp = bsxfun(@times,kSpace_update_temp,mask);

    for i=1:nSMS
        kSpace_update_radial(:,:,:,:,:,:,i) = GROG.GROG_cart2rad_4P(kSpace_update_temp(:,:,:,:,:,:,i),G{i});
        fidelity_update_radial(:,:,:,:,:,:,i) = kSpace_radial(:,:,:,:,:,:,i) - kSpace_update_radial(:,:,:,:,:,:,i); 
        fidelity_update(:,:,:,:,:,:,i) = GROG.GROG_rad2cart_4P(fidelity_update_radial(:,:,:,:,:,:,i),G{i});
    end

    
    fidelity_update_temp(:,:,:,:,1) = sum(fidelity_update,7);
    fidelity_update_temp(:,:,:,:,2) = fidelity_update(:,:,:,:,:,:,1) + exp(1i*2*pi/3)*fidelity_update(:,:,:,:,:,:,2) + exp(1i*4*pi/3)*fidelity_update(:,:,:,:,:,:,3);
    fidelity_update_temp(:,:,:,:,3) = fidelity_update(:,:,:,:,:,:,1) + exp(1i*4*pi/3)*fidelity_update(:,:,:,:,:,:,2) + exp(1i*2*pi/3)*fidelity_update(:,:,:,:,:,:,3);
    %fidelity_update_temp = fftshift(fidelity_update_temp,1);
    %fidelity_update_temp = fftshift(fidelity_update_temp,2);
    fidelity_update_temp = bsxfun(@times,fidelity_update_temp,fftshift_mask);
    fidelity_update = ifft2(fidelity_update_temp); clear 
    %fidelity_update = bsxfun(@times,fidelity_update,fftshift_mask);
    %fidelity_update = fftshift(fidelity_update,1);
    %fidelity_update = fftshift(fidelity_update,2);
%}    
    clear kSpace_update kSpace_update_radial fidelity_update_radial
    fidelity_update = sum(bsxfun(@times,fidelity_update,sens_map_conj),4);
    
    para.CPUtime.fidelity(iter_no) = toc;
    tic; tTV_update = compute_tTV_yt(new_img_x,weight_tTV,beta_sqrd);  para.CPUtime.tTV(iter_no) = toc;
    tic; sTV_update = compute_sTV_yt(new_img_x,weight_sTV,beta_sqrd);  para.CPUtime.sTV(iter_no) = toc;
    %tic,sTV_update = weight_sTV*TVG_nD_matrix(new_img_x,beta_sqrd);para.CPUtime.sTV(iter_no) = toc;
    update_term  = fidelity_update + tTV_update + sTV_update; clear tTV_update sTV_update
    
    tic;
    if debug
        [fidelityNorm(iter_no), temporalNorm(iter_no), spatialNorm(iter_no), totalCost(iter_no)] = NormCalculation(fidelity_update, new_img_x, weight_sTV, weight_tTV);
        if iter_no>1 && ifBLS==1
            if totalCost(iter_no)>totalCost(iter_no-1)
                step_size = step_size*0.5;
                %disp('Chaning step size into');disp(step_size);
            end
        end 
    end
    clear fidelity_update
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

    %if ifplot
        showImage(new_img_x,fidelityNorm,temporalNorm,spatialNorm,totalCost)
    %end
    
    if mod(iter_no,save_frequency)==0
        Image = squeeze(new_img_x);
        disp('Saving image into Results...')
        save_dir_temp = strcat(para.Recon.save_dir,para.time);
        save([save_dir_temp,'.mat'],'Image','new_img_z','para','-v7.3');
        if debug
            save([save_dir_temp,'.mat'],'temporalNorm','spatialNorm','fidelityNorm','totalCost','-append');
        end
    end
    
    toc(t1);
end

Image = squeeze(new_img_x);
para = get_CPU_time(para);
fprintf(['Iterative STCR running time is ' num2str(para.CPUtime.interative_recon) 's' '\n'])
disp('Reconstruction done');fprintf('\n')
disp('Saving image into Results...')
Image = gather(Image);
sx = size(Image,1);
fftshift_mask = single(ones(sx));
fftshift_mask(1:2:sx,:) = -fftshift_mask(1:2:sx,:);
fftshift_mask(:,1:2:sx) = -fftshift_mask(:,1:2:sx);
Image = bsxfun(@times,Image,-fftshift_mask);
save_dir = strcat(para.Recon.save_dir,para.time);
save([save_dir,'.mat'],'Image','para','-v7.3');
if debug
    temporalNorm = gather(temporalNorm);
    spatialNorm = gather(spatialNorm);
    fidelityNorm = gather(fidelityNorm);
    totalCost = gather(totalCost);
    save([save_dir,'.mat'],'temporalNorm','spatialNorm','fidelityNorm','totalCost','-append')
end
disp('All done.');fprintf('\n')