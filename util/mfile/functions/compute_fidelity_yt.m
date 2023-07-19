function fidelity_update = compute_fidelity_yt(ifNUFFT,nSMS,kSpace_type,image,kSpace,sens_map,sens_map_conj,phase_mod,phase_mod_conj,varargin)

if ifNUFFT == 1

    NUFFT = varargin{1}{1};
    %NUFFT_FE = varargin{2};
    kSpace_update = NUFFT * bsxfun(@times,image,sens_map);
    %kSpace_update = NUFFT_FE * bsxfun(@times,image,sens_map);
    if nSMS ~= 1
        kSpace_update = kSpace - sum(bsxfun(@times,kSpace_update,phase_mod),5);
        kSpace_update = bsxfun(@times,kSpace_update,phase_mod_conj);
    else
        kSpace_update = kSpace - kSpace_update;
    end

    fidelity_update = NUFFT'* kSpace_update;
    %fidelity_update = NUFFT'* kSpace_update(:,RayIndex_A3J,:,:,:,:);
    %fidelity_update = fidelity_update*7.9134; % number need to test
    %fidelity_update = kSpace - fidelity_update;
        
else

    kSpace_update = bsxfun(@times,image,sens_map);
    kSpace_update(:,:,:,1:2,:,:) = fft2(kSpace_update(:,:,:,1:2,:,:)); % I think only 1% difference in time if split into 2 on CPU.
    kSpace_update(:,:,:,3:4,:,:) = fft2(kSpace_update(:,:,:,3:4,:,:));
    kSpace_update(:,:,:,5:6,:,:) = fft2(kSpace_update(:,:,:,5:6,:,:));
    kSpace_update(:,:,:,7:end,:,:) = fft2(kSpace_update(:,:,:,7:end,:,:));
        
    if nSMS~=1
        switch kSpace_type
            case 'aliasing_est'
                fidelity_update = bsxfun(@times,kSpace_update,phase_mod(:,:,:,1,1)); clear kSpace_update;
            otherwise
                kSpace_update = sum(bsxfun(@times,kSpace_update,phase_mod),5);
        end
    end
    
    switch kSpace_type
        case 'imspace'
            fidelity_update = bsxfun(@times,kSpace_update,phase_mod_conj); clear kSpace_update;
            fidelity_update(:,:,:,1:2,:,:) = ifft2(fidelity_update(:,:,:,1:2,:,:));
            fidelity_update(:,:,:,3:4,:,:) = ifft2(fidelity_update(:,:,:,3:4,:,:));
            fidelity_update(:,:,:,5:6,:,:) = ifft2(fidelity_update(:,:,:,5:6,:,:));
            fidelity_update(:,:,:,7:8,:,:) = ifft2(fidelity_update(:,:,:,7:8,:,:));
            
            fidelity_update = kSpace - fidelity_update;
        case 'kspace'
            fidelity_update = kSpace - kSpace_update; clear kSpace_update;
            fidelity_update = bsxfun(@times,fidelity_update,phase_mod_conj);
            fidelity_update(:,:,:,1:2,:,:) = ifft2(fidelity_update(:,:,:,1:2,:,:));
            fidelity_update(:,:,:,3:4,:,:) = ifft2(fidelity_update(:,:,:,3:4,:,:));
            fidelity_update(:,:,:,5:6,:,:) = ifft2(fidelity_update(:,:,:,5:6,:,:));
            fidelity_update(:,:,:,7:end,:,:) = ifft2(fidelity_update(:,:,:,7:end,:,:));
        case 'aliasing_est'
            fidelity_update(:,:,:,1:2,:,:) = ifft2(fidelity_update(:,:,:,1:2,:,:));
            fidelity_update(:,:,:,3:4,:,:) = ifft2(fidelity_update(:,:,:,3:4,:,:));
            fidelity_update(:,:,:,5:6,:,:) = ifft2(fidelity_update(:,:,:,5:6,:,:));
            fidelity_update(:,:,:,7:8,:,:) = ifft2(fidelity_update(:,:,:,7:8,:,:));
            
            fidelity_update = kSpace - fidelity_update;
    end
end

fidelity_update = sum(bsxfun(@times,fidelity_update,sens_map_conj),4);
    
end