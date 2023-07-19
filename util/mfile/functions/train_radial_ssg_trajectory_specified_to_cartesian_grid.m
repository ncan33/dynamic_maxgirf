function kernel = train_radial_ssg_trajectory_specified_to_cartesian_grid(im_sb, theta, kx, ky, phase_mod, k_test)

fprintf([repmat('-', [1, 75]), '\n'])
fprintf(sprintf('begin trajectory specified radial slice GRAPPA calibration: \n'));
fprintf([repmat('-', [1, 75]), '\n'])
                                                                           
tic

kspace_segments = 4;
kernel_length = 5;

point_extend = (kernel_length - 1)/2;

[sx, ~, nof_calib, nc, nsms] = size(im_sb);
nor_target = size(kx, 2);
nof_target = size(kx, 3);


phase_mod = reshape(phase_mod, [1, nor_target, nof_target]);
phase_mod_conj = exp(-1i*2*pi*phase_mod/3);
phase_mod_conj = cat(5, ones(size(phase_mod_conj)), phase_mod_conj, conj(phase_mod_conj));

theta = reshape(theta, [1, nor_target, nof_target]);

kernel = zeros([nc*nsms, nc*nsms*kernel_length, kspace_segments, nor_target, nof_target], 'like', im_sb);

% get which cartesian k-space to estimate
kx_cart = round(kx + sx/2+1);
ky_cart = round(ky + sx/2+1);
kx_cart(kx_cart < 1) = sx - kx_cart(kx_cart < 1);
ky_cart(ky_cart < 1) = sx - ky_cart(ky_cart < 1);
kx_cart(kx_cart > sx) = kx_cart(kx_cart > sx) - sx;
ky_cart(ky_cart > sx) = ky_cart(ky_cart > sx) - sx;

% get Cartesian k-space
kSpace_cart_sb = fftshift2(fft2(fftshift2(im_sb)));

% display some text
text_display = sprintf('number of time frames: %16.0f\nnumber of calibration time frames: %4.0f\nnumber of coils: %22.0f\nnumber of slices: %21.0f\n\n', nof_target, nof_calib, nc, nsms);
fprintf(text_display)
text_display = sprintf('kernel length: %24.0f\nnumber of k-space segments: %11.0f\n\n', kernel_length, kspace_segments);
fprintf(text_display)

for iframe = 1:nof_target
if iframe ~= 1
    fprintf(repmat('\b',1,linelength));
end
linelength = fprintf(sprintf('calibrate frame %g / %g \n', iframe, nof_target));

% get trajectory specified calibration data
% interpolate single band image into target trajectory
kx_temp = kx(:, :, iframe);
ky_temp = ky(:, :, iframe);
N = NUFFT.init(kx_temp, ky_temp, 1, [6, 6], sx, sx);

kSpace_sb = zeros([sx, nor_target, nof_calib, nc, nsms], 'like', im_sb);
for itemp = 1:nof_calib
    kSpace_sb(:, :, itemp, :, :) = NUFFT.NUFFT(im_sb(:, :, itemp, :, :), N);
end

idx_cart = sub2ind([sx, sx], kx_cart(:, :, iframe), ky_cart(:, :, iframe));
idx_cart = repmat(idx_cart, [1, 1, nof_calib, nc, nsms]);
idx_cart = idx_cart + reshape((0:nof_calib-1)*sx*sx, [1, 1, nof_calib]);
idx_cart = idx_cart + reshape((0:nc-1)*sx*sx*nof_calib, [1, 1, 1, nc]);
idx_cart = idx_cart + reshape((0:nsms-1)*sx*sx*nof_calib*nc, [1, 1, 1, 1, nsms]);
kSpace_cart_target = kSpace_cart_sb(idx_cart);

phase_mod_temp = conj(phase_mod_conj(:, :, iframe, :, :));
phase_idx_temp = phase_mod(:, :, iframe);
kSpace_mb = sum(kSpace_sb .* phase_mod_temp, 5);
theta_temp = theta(:, :, iframe);

for iray = 1:nor_target
    iray
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
    
    kSpace_source_phase_mod_1 = kSpace_mb(:, ray_idx_phase_mod_1, :, :);
    kSpace_source_phase_mod_2 = kSpace_mb(:, ray_idx_phase_mod_2, :, :);
    
    % this doesn't work for 2*pi 
    % flip the ray under some condition.
    % When the rays are close around 0-2*pi, or pi, the rays are actually
    % close but you need to flip them
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
    
    
    source = cat(2, kSpace_mb(:, iray, :, :, :), kSpace_source_phase_mod_1, kSpace_source_phase_mod_2);
    target = kSpace_sb(:, iray, :, :, :);
    %{
    for iseg = 1:sx
        if iseg <= point_extend
            segment_idx = 1:kernel_length;
        elseif iseg > sx - point_extend
            segment_idx = sx - kernel_length + 1 : sx;
        else
            segment_idx = iseg - point_extend : iseg + point_extend;
        end
        
        source_temp = source(segment_idx, :, :, :);
        target_temp = target(iseg, :, :, :, :);
        
        % the nsms here is also number of calibration rays
        source_temp = permute(source_temp, [1, 2, 4, 3]);
        source_temp = reshape(source_temp, [kernel_length*nsms*nc, nof_calib]);
        target_temp = permute(target_temp, [4, 5, 3, 1, 2]);
        target_temp = reshape(target_temp, [nc*nsms, nof_calib]);
        
        % add noise into calibration
%         noise_level = max(abs(source_kernel(:)))/50;
%         source_noise = randn(size(source_kernel), 'single') + 1i * randn(size(source_kernel), 'single');
%         source_noise = noise_level * source_noise;
%         source_kernel = cat(2, source_kernel, source_noise);
%         target_kernel = cat(2, target_kernel, zeros(size(target_kernel), 'single'));
        % end add noise
        kernel(:, :, iseg, iray, iframe) = target_temp * pinv(source_temp);
%         k_test(:, :, iseg, iray) = kernel(:, :, iseg, iray, iframe) * source_temp;
    end
    %}
    for iseg = 1:kspace_segments
        segment_idx = (1:sx/kspace_segments) + (iseg-1)*sx/kspace_segments;
        if iseg ~= 1
            segment_idx = [segment_idx(1) - (1:point_extend), segment_idx];
        end
        if iseg ~= kspace_segments
            segment_idx = [segment_idx, segment_idx(end) + (1:point_extend)];
        end
        
        source_temp = source(segment_idx, :, :, :, :);
        target_temp = target(segment_idx, :, :, :);
        nkernel = (length(segment_idx) - kernel_length + 1) * nof_calib;
        
        % the nsms here is also number of calibration rays
        source_kernel = zeros([kernel_length * nc * nsms, nkernel], 'single');
        target_kernel = zeros([nc*nsms, nkernel], 'single');
        
        for itime = 1:nof_calib
            for ipatch = 1:nkernel/nof_calib
                point_idx = (1:kernel_length) + (ipatch-1);
                center = point_idx((kernel_length + 1)/2);
                source_kernel_temp = source_temp(point_idx, :, itime, :);
                source_kernel(:, (itime-1)*nkernel/nof_calib + ipatch) = source_kernel_temp(:);
                target_kernel(:, (itime-1)*nkernel/nof_calib + ipatch) = target_temp(center, :, itime, :);
            end
        end
        
        % add noise into calibration
%         noise_level = max(abs(source_kernel(:)))/50;
%         source_noise = randn(size(source_kernel), 'single') + 1i * randn(size(source_kernel), 'single');
%         source_noise = noise_level * source_noise;
%         source_kernel = cat(2, source_kernel, source_noise);
%         target_kernel = cat(2, target_kernel, zeros(size(target_kernel), 'single'));
        % end add noise
        kernel(:, :, iseg, iray, iframe) = target_kernel * pinv(source_kernel);
    end
    
end

end

% k_test = reshape(k_test, [nc, nsms, nof_calib, sx, nor_target]);
% k_test = permute(k_test, [4, 5, 3, 1, 2]);
% k_test_cart = zeros([sx, sx, nof_calib, nc, nsms]);
% k_test_cart(idx_cart) = k_test;
% im_test = fftshift2(ifft2(fftshift2(k_test_cart)));

t1 = toc;
fprintf(repmat('\b',1,linelength));
fprintf(sprintf('calibration done, time = %.2f s\n', t1));
fprintf([repmat('-', [1, 75]), '\n'])

k_temp = apply_radial_ssg_trajectroy_speficied_to_cartesian_grid(k_test, phase_mod, theta, kernel);
k_test_cart = zeros([sx, sx, 1, nc, nsms]);
k_test_cart(idx_cart(:,:,1,:,:)) = k_temp;

keyboard