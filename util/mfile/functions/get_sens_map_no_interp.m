function [sens_map,para] = get_sens_map_no_interp(kSpace,phase_mod,para)

sens_map_dir = para.Recon.sens_map_dir;
noSensMap    = para.Recon.noSensMap;

if noSensMap == 0
    disp('Senstivity map already exsists...');tic
    load([pwd,sens_map_dir])
else
    disp('Estimating senstivity map...');tic

    nSMS          = para.Recon.nSMS;
    [sx,~,~,no_comp,~,ns] = size(kSpace); 
    
        if nSMS>1
            img_forSensMap = ifft2(bsxfun(@times,kSpace,conj(phase_mod)));
        else
            img_forSensMap = ifft2(kSpace);
        end
        
        img_forSensMap = squeeze(sum(img_forSensMap,3));

    img_forSensMap = reshape(img_forSensMap,[sx sx no_comp nSMS*ns]);
    sens_map = single(zeros(sx,sx,1,no_comp,nSMS*ns));
    for i=1:nSMS*ns
        sens_map(:,:,1,:,i) = ismrm_estimate_csm_walsh_optimized_yt(img_forSensMap(:,:,:,i),10);
        %sens_map(:,:,1,:,i) = ismrm_estimate_csm_mckenzie(img_forSensMap(:,:,:,i));
        %sens_map(:,:,1,:,i) = ir_mri_sensemap_admm(img_forSensMap(:,:,:,i));
        %[~,sens_map(:,:,1,:,i)] = adapt_array_2d(img_forSensMap(:,:,:,i));
        % block wise, maybe not so good?
    end
    sens_map = reshape(sens_map,[sx sx 1 no_comp nSMS ns]);
    sens_map_scale = max(abs(sens_map(:)));                                
    sens_map = sens_map/sens_map_scale;
    sens_map_conj = conj(sens_map);
    sens_correct_term = 1./sum(sens_map_conj.*sens_map,4);
    sens_correct_term = sqrt(sens_correct_term);
    sens_map = bsxfun(@times,sens_correct_term,sens_map);
    %save([pwd,sens_map_dir],'sens_map');
end

%if isfield(para,'slice_pick')
%    sens_map = sens_map(:,:,:,:,:,slice_pick);
%end

para.CPUtime.estimate_sens_map = toc;toc;fprintf('\n');
