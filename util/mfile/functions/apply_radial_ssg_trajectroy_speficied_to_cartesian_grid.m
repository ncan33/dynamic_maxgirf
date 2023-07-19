function kSpace_sb = apply_radial_ssg_trajectroy_speficied_to_cartesian_grid(kSpace_mb, phase_mod, theta, kernel)

[sx, nor, nof, nc] = size(kSpace_mb);
phase_mod = reshape(phase_mod, [1, nor, nof]);
theta = reshape(theta, [1, nor, nof]);
[~, nsource, kspace_segments, nor, nof] = size(kernel);
nsms = 3;
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
    
    
    for iseg = 1:sx
        if iseg <= point_extend
            segment_idx = 1:kernel_length;
        elseif iseg > sx - point_extend
            segment_idx = sx - kernel_length + 1 : sx;
        else
            segment_idx = iseg - point_extend : iseg + point_extend;
        end
        
        source_temp = source(segment_idx, :, :, :);
        % the nsms here is also number of calibration rays
        source_temp = permute(source_temp, [1, 2, 4, 3]);
        source_temp = reshape(source_temp, [kernel_length*nsms*nc, nof]);

        kspace_temp = kernel(:, :, iseg, iray, iframe) * source_temp;
        kspace_temp = reshape(kspace_temp, [nc, nsms]);
        kSpace_sb(iseg, iray, iframe, :, :) = kspace_temp;
    end
end
end
end
% fprintf(repmat('\b',1,linelength));