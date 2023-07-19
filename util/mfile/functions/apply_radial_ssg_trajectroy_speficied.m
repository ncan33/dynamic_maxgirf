function kSpace_sb = apply_radial_ssg_trajectroy_speficied(kSpace_mb, phase_mod, theta, kernel)

[sx, nor, nof, nc] = size(kSpace_mb);
phase_mod = reshape(phase_mod, [1, nor, nof]);
theta = reshape(theta, [1, nor, nof]);
[~, nsource, kspace_segments, nor, nsms, nof] = size(kernel);
kernel_length = nsource/nc/nsms;
point_extend = (kernel_length - 1)/2;

kSpace_sb = zeros([sx, nor, nof, nc, nsms], 'like', kSpace_mb);
% test_frame = 1:100;
% ntest_frame = length(test_frame);
% kSpace_temp = zeros([sx, nor_target, ntest_frame, nc, islice]);
% test_kspace = reshape(test_kspace, [sx, nor_target, nof_target, nc]);
% kSpace_mb = test_kspace(:, :, test_frame, :);

for iframe = 1:nof
% if iframe ~= 1
%     fprintf(repmat('\b',1,linelength));
% end
% linelength = fprintf(sprintf('apply kernel to frame %g / %g \n', iframe, nof));

phase_idx_temp = phase_mod(:, :, iframe);
theta_temp = theta(:, :, iframe);

for islice = 1:nsms
for iray = 1:nor
    phase_idx = phase_idx_temp(1, iray);
    if phase_idx == 0
        phase_idx_1 = 1;
        phase_idx_2 = 2;
    elseif phase_idx == 1
        phase_idx_1 = 0;
        phase_idx_2 = 2;
    elseif phase_idx == 2
        phase_idx_1 = 0;
        phase_idx_2 = 1;
    end
    
    theta_this = theta_temp(1, iray);
    
    phase_mod_1_idx = find(phase_idx_temp == phase_idx_1);
    phase_mod_2_idx = find(phase_idx_temp == phase_idx_2);
    
    theta_phase_mod_1 = theta_temp(:, phase_mod_1_idx);
    theta_phase_mod_2 = theta_temp(:, phase_mod_2_idx);
    
    dtheta_phase_mod_1 = mod(theta_this - theta_phase_mod_1, pi);
    dtheta_phase_mod_1 = min(dtheta_phase_mod_1, pi - dtheta_phase_mod_1);
    dtheta_phase_mod_2 = mod(theta_this - theta_phase_mod_2, pi);
    dtheta_phase_mod_2 = min(dtheta_phase_mod_2, pi - dtheta_phase_mod_2);
    
    [min_dtheta_phase_mod_1, ray_idx_phase_mod_1] = min(dtheta_phase_mod_1);
    ray_idx_phase_mod_1 = phase_mod_1_idx(ray_idx_phase_mod_1);
    
    [min_dtheta_phase_mod_2, ray_idx_phase_mod_2] = min(dtheta_phase_mod_2);
    ray_idx_phase_mod_2 = phase_mod_2_idx(ray_idx_phase_mod_2);
    
    kSpace_source_phase_mod_1 = kSpace_mb(:, ray_idx_phase_mod_1, iframe, :);
    kSpace_source_phase_mod_2 = kSpace_mb(:, ray_idx_phase_mod_2, iframe, :);
    
    % this doesn't work for 2*pi 
    if theta_this < min_dtheta_phase_mod_1 && pi - theta_temp(ray_idx_phase_mod_1) < min_dtheta_phase_mod_1
        kSpace_source_phase_mod_1 = flipud(kSpace_source_phase_mod_1);
    end
    if theta_this < min_dtheta_phase_mod_2 && pi - theta_temp(ray_idx_phase_mod_2) < min_dtheta_phase_mod_2
        kSpace_source_phase_mod_2 = flipud(kSpace_source_phase_mod_2);
    end
    if pi - theta_this < min_dtheta_phase_mod_1 && theta_temp(ray_idx_phase_mod_1) < min_dtheta_phase_mod_1
        kSpace_source_phase_mod_1 = flipud(kSpace_source_phase_mod_1);
    end
    if pi - theta_this < min_dtheta_phase_mod_2 && theta_temp(ray_idx_phase_mod_2) < min_dtheta_phase_mod_2
        kSpace_source_phase_mod_2 = flipud(kSpace_source_phase_mod_2);
    end
    
    
    source = cat(2, kSpace_mb(:, iray, iframe, :), kSpace_source_phase_mod_1, kSpace_source_phase_mod_2);
%     target = kSpace_sb(:, i, :, :, islice);
    
    
    for iseg = 1:kspace_segments
        segment_idx = (1:sx/kspace_segments) + (iseg-1)*sx/kspace_segments;
        if iseg ~= 1
            segment_idx = [segment_idx(1) - (1:point_extend), segment_idx];
        end
        if iseg ~= kspace_segments
            segment_idx = [segment_idx, segment_idx(end) + (1:point_extend)];
        end
        source_temp = source(segment_idx, :, :, :);
%         target_temp = target(segment_idx, :, :, :);
        nkernel = (length(segment_idx) - kernel_length + 1);
        
        % the nsms here is also number of calibration rays
        source_kernel = zeros([kernel_length * nc * nsms, nkernel], 'like', kSpace_mb);
%         target_kernel = zeros([nc, nkernel], 'single');
        
        for ipatch = 1:nkernel
            point_idx = (1:kernel_length) + (ipatch-1);
%             center = point_idx((kernel_length + 1)/2);
            source_kernel_temp = source_temp(point_idx, :, :, :);
            source_kernel(:, ipatch) = source_kernel_temp(:);
%             target_kernel(:, (itime-1)*nkernel/nof + ipatch) = target_temp(center, :, itime, :);
        end

        kspace_temp = kernel(:, :, iseg, iray, islice, iframe) * source_kernel;
        kSpace_sb(segment_idx(point_extend+1:end-point_extend), iray, iframe, :, islice) = kspace_temp.';
    end
end
end
end
% fprintf(repmat('\b',1,linelength));