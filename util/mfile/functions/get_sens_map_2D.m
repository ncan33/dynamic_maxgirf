function [sens_map,para] = get_sens_map(kSpace_all,theta_all,phase_mod_all,RayPosition,para)

sens_map_dir = para.Recon.sens_map_dir;
noSensMap    = para.Recon.noSensMap;

if noSensMap == 0
    disp('Senstivity map already exsists...');tic
    load([pwd,sens_map_dir])
else
    disp('Estimating senstivity map...');tic
    
    nor_SM        = para.Recon.nor_SM;
    nof           = para.Recon.nof;
    nSMS          = para.Recon.nSMS;
    kCenter       = para.kSpace_center;
    ifNUFFT       = para.ifNUFFT;
    ifasym        = para.asymmetry;
    
    RayIndex_SM   = RayPosition < 1 + nor_SM;
    
    [sx,~,no_comp,~,ns] = size(kSpace_all);
    kSpace_SM = kSpace_all(:,RayIndex_SM,:,:,:);
    kSpace_SM = reshape(kSpace_SM,[sx,nor_SM,nof,no_comp,1,ns]);
    theta_SM = theta_all(RayIndex_SM);
    theta_SM = reshape(theta_SM,[1,nor_SM,nof]);    
    phase_mod_SM = phase_mod_all(1,RayIndex_SM,1,:);
    phase_mod_SM = reshape(phase_mod_SM,[1,nor_SM,nof,1,nSMS]);

    if ifNUFFT ==1
        kSpace_SM = bsxfun(@times,kSpace_all,phase_mod_all);
        [x_coor_SM, y_coor_SM] = get_k_coor(sx,theta_all,ifNUFFT,kCenter);
        coor_SM = x_coor_SM + 1i*y_coor_SM;
        if ifasym
            coor_SM(asymRayStart:asymRayEnd,:,:) = [];
            kSpace_SM(asymRayStart:asymRayEnd,:,:,:,:) = [];
        end
        center_shift = kCenter - sx/2 -1;
        W = designFilter(sx,center_shift,'hamming');
        NUFFT_SM = FftTools.MultiNufft(coor_SM,W,0,[sx sx]);
        NUFFT_SM.kSize_z = no_comp;
        img_forSensMap = NUFFT_SM'* kSpace_SM;
    else
        interp_method = para.interp_method;
        
        [x_coor_SM,y_coor_SM] = get_k_coor(sx,theta_SM,ifNUFFT,kCenter);
        switch interp_method
            case 'grid3'
                kSpace_SM = pre_interp_radial_to_cart(kSpace_SM,x_coor_SM,y_coor_SM);
            case 'nn'
                kSpace_SM = permute(kSpace_SM,[1 2 4 3 5 6]);
                kSpace_c = zeros([sx sx no_comp nof 1 ns]);
                parfor i=1:ns
                    kSpace_c(:,:,:,:,1,i) = interp_nn(kSpace_SM(:,:,:,:,1,i), x_coor_SM,y_coor_SM);
                end
                kSpace_SM = permute(kSpace_c,[1 2 4 3 5 6]); clear kSpace_c
            case 'GROG'
                kSpace = GROG.GROG(squeeze(kSpace),x_coor,y_coor, 0, 1);
                kSpace = permute(kSpace,[1 2 3 4 6 5]);
                %kSpace_SM = permute(kSpace_SM,[1 2 4 3 5 6]);
                %kSpace_c = zeros([sx sx no_comp nof 1 ns]);
                %parfor i=1:ns*nSMS
                %    kSpace_c(:,:,:,:,1,i) = scGROG(kSpace_SM(:,:,:,:,1,i), x_coor_SM,y_coor_SM);
                %end
                %kSpace_SM = permute(kSpace_c,[1 2 4 3 5 6]); clear kSpace_c
        end
        
        fftshift_mask = single(ones(sx));
        fftshift_mask(1:2:sx,:) = -fftshift_mask(1:2:sx,:);
        fftshift_mask(:,1:2:sx) = -fftshift_mask(:,1:2:sx);
        kSpace_SM = bsxfun(@times,kSpace_SM,fftshift_mask);
        
        if nSMS>1
            phase_mod_SM_all = repmat(phase_mod_SM,[sx,1,1,1,1]); clear phase_mod_SM
            phase_mod_SM = single(logical(abs(kSpace_SM(:,:,:,1,1,1))));
            phase_mod_SM(:,:,:,:,2:nSMS) = pre_interp_radial_to_cart(phase_mod_SM_all(:,:,:,:,2:nSMS),x_coor_SM,y_coor_SM);
            phase_mod_SM_conj = conj(phase_mod_SM);
            img_forSensMap = ifft2(bsxfun(@times,kSpace_SM,phase_mod_SM_conj));
        else
            img_forSensMap = ifft2(kSpace_SM);
        end
        
        img_forSensMap = squeeze(sum(img_forSensMap,3));
    end

    img_forSensMap = reshape(img_forSensMap,[sx sx no_comp nSMS*ns]);
    sens_map = single(zeros(sx,sx,1,no_comp,nSMS*ns));
    parfor i=1:nSMS*ns
        sens_map(:,:,1,:,i) = ismrm_estimate_csm_walsh_optimized_yt(img_forSensMap(:,:,:,i));
    end
    sens_map = reshape(sens_map,[sx sx 1 no_comp nSMS ns]);
    sens_map_scale = max(abs(sens_map(:)));                                
    sens_map = sens_map/sens_map_scale;
    sens_map_conj = conj(sens_map);
    sens_correct_term = 1./sum(sens_map_conj.*sens_map,4);
    sens_correct_term = sqrt(sens_correct_term);
    sens_map = bsxfun(@times,sens_correct_term,sens_map);
    save([pwd,sens_map_dir],'sens_map');
end

if isfield(para,'slice_pick')
    sens_map = sens_map(:,:,:,:,:,para.slice_pick);
end

para.CPUtime.estimate_sens_map = toc;toc;fprintf('\n');
