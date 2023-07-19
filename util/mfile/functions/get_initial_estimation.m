function [first_est,para] = get_initial_estimation(sens_map,kSpace_all,theta_all,phase_mod_all,RayPosition,kSpace,phase_mod,para)
disp('Calculating initial estimation...');tic

if para.FirstEstimationA3J

    selected_rays_start_FE = para.selected_rays_start_FE;
    selected_rays_end_FE   = para.selected_rays_end_FE;
    RayIndex_FE = RayPosition > selected_rays_start_FE-1 .* RayPosition < selected_rays_end_FE+1;
    nor_FE = selected_rays_end_FE - selected_rays_start_FE + 1;
    
    theta_A3J = theta_all(RayIndex_FE);
    theta_A3J = reshape(theta_A3J,[1,nor,nof]);
    theta_A3J_temp = zeros([1,nor_FE*3,nof]);
    theta_A3J_temp(:,nor_FE+1:2*nor_FE,:) = theta_A3J;
    theta_A3J_temp(:,1:nor_FE,2:end) = theta_A3J(:,:,1:end-1);
    theta_A3J_temp(:,2*nor_FE+1:end,1:end-1) = theta_A3J(:,:,2:end);
    theta_A3J = theta_A3J_temp; clear theta_A3J_temp
    
    phase_mod_FE = phase_mod_all(1,RayIndex_FE,1,:);
    phase_mod_FE = reshape(phase_mod_FE,[1,nor,nof,1,nSMS]);
    phase_mod_temp = single(zeros([1,nor_FE*3,nof,1,3]));
    phase_mod_temp(:,nor_FE+1:2*nor_FE,:,:,:) = phase_mod_FE;
    phase_mod_temp(:,1:nor_FE,2:end,:,:) = phase_mod_FE(:,:,1:end-1,:,:);
    phase_mod_temp(:,2*nor_FE+1:end,1:end-1,:,:) = phase_mod_FE(:,:,2:end,:,:);
    phase_mod_FE = phase_mod_temp; clear phase_mod2_temp
    
    kSpace_FE = kSpace_all(:,RayIndex_FE,:,:,:);
    kSpace_FE = reshape(kSpace_FE,[sx,nor,nof,no_comp,1,ns]);
    kSpace_FE_temp = single(zeros([sx nor_FE*3 nof no_comp 1 ns]));
    kSpace_FE_temp(:,nor_FE+1:2*nor_FE,:,:,:,:) = kSpace_FE;
    kSpace_FE_temp(:,1:nor_FE,2:end,:,:,:) = kSpace_FE(:,:,1:end-1,:,:,:);
    kSpace_FE_temp(:,2*nor_FE+1:end,1:end-1,:,:,:) = kSpace_FE(:,:,2:end,:,:,:);
    kSpace_FE = kSpace_FE_temp; clear kSpace_FE_temp 

        [x_coor_FE,y_coor_FE] = get_k_coor(sx,theta_A3J,ifNUFFT,kCenter);
        kSpace_FE = pre_interp_radial_to_cart(kSpace_FE,x_coor_FE,y_coor_FE);
        kSpace_FE = bsxfun(@times,kSpace_FE,fftshift_mask);
        phase_mod_FE_all = repmat(phase_mod_FE,[sx,1,1,1,1]); clear phase_mod_FE
        phase_mod_FE(:,:,:,:,1)   = pre_interp_radial_to_cart(phase_mod_FE_all(:,:,:,:,1),x_coor_FE,y_coor_FE);
        phase_mod_FE(:,:,:,:,2:3) = pre_interp_radial_to_cart(phase_mod_FE_all(:,:,:,:,2:3),x_coor_FE,y_coor_FE);
        phase_mod_FE = fftshift(fftshift(phase_mod_FE,1),2);
        phase_mod_FE_conj = conj(phase_mod_FE);
        img_FE = ifft2(bsxfun(@times,kSpace_FE,phase_mod_FE_conj));
        %img_k = ifft2(bsxfun(@times,kSpace,phase_mod_conj));
else

    img_FE = ifft2(bsxfun(@times,kSpace,conj(phase_mod)));

    %img_k = ifft2(kSpace);
    %img_k = img_FE;
end

img_FE = img_FE/mean(abs(img_FE(:)))*0.0215;
%img_k = img_k/mean(abs(img_k(:)))*0.0215;

first_est = single(sum(bsxfun(@times, img_FE, conj(sens_map)),4));

%first_est = repmat(sum(first_est,3),[1 1 nof 1 1 1])/nof;

para.CPUtime.initial_est = toc;toc;fprintf('\n');