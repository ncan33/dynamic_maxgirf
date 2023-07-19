function [sens_map,para] = get_sens_map_3D(kSpace,para)

%sens_map_dir = para.Recon.sens_map_dir;
%noSensMap    = para.Recon.noSensMap;
kSpace_dim   = length(size(kSpace));

%if noSensMap == 0
%    disp('Senstivity map already exsists...');tic
%    load([pwd,sens_map_dir])
%else
    disp('Estimating senstivity map...');tic
    if length(size(kSpace)) == 4
        [sx,sy,sz,no_comp] = size(kSpace);
        img_forSensMap = ifft3(kSpace);
    else
        [sx,sy,sz,~,no_comp] = size(kSpace);
        img_forSensMap = squeeze(sum(ifft3(kSpace),4));
    end
%     [sx,sy,sz,~,~] = size(kSpace);
%fftshift_mask = single(ones([sx sy sz]));
%fftshift_mask(1:2:sx,:,:) = -fftshift_mask(1:2:sx,:,:);
%fftshift_mask(:,1:2:sy,:) = -fftshift_mask(:,1:2:sy,:);
%fftshift_mask(:,:,1:2:sz) = -fftshift_mask(:,:,1:2:sz);
%    for i=1:52; for j=1:26; img_forSensMap0(:,:,i,j) = imresize(img_forSensMap(:,:,i,j).*fftshift_mask(:,:,1,1),[640 320]); end; end
%    img_forSensMap = img_forSensMap0; clear img_forSensMap0; sy = 320;
%    keyboard
    sens_map = single(zeros(sx,sy,sz,1,no_comp));
    for i=1:sz
        sens_map(:,:,i,1,:) = ismrm_estimate_csm_walsh_optimized_yt(squeeze(img_forSensMap(:,:,i,:)),10);
    end
    if length(size(kSpace)) == 4
        sens_map = squeeze(sens_map);
    end
    sens_map_scale = max(abs(sens_map(:)));                                
    sens_map = sens_map/sens_map_scale;
    sens_map_conj = conj(sens_map);
    sens_correct_term = 1./sum(sens_map_conj.*sens_map,kSpace_dim);
    sens_correct_term = sqrt(sens_correct_term);
    sens_map = bsxfun(@times,sens_correct_term,sens_map);
    %save([pwd,sens_map_dir],'sens_map','-v7.3');
%end

para.Recon.kSpace_dim = kSpace_dim;
para.CPUtime.estimate_sens_map = toc;toc;fprintf('\n');
