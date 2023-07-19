function kernel = train_radial_ssg_trajectory_specified(im_sb, theta, kx, ky, phase_mod, sens, test_kspace)

fprintf([repmat('-', [1, 75]), '\n'])
fprintf(sprintf('begin trajectory specified radial slice GRAPPA calibration: \n'));
fprintf([repmat('-', [1, 75]), '\n'])
                                                                           
tic
% nor_target = 30;
kspace_segments = 4;
kernel_length = 5;

point_extend = (kernel_length - 1)/2;

[sx, ~, nof_calib, nc, nsms] = size(im_sb);
nor_target = size(kx, 2);
nof_target = size(kx, 3);


phase_mod = reshape(phase_mod, [1, nor_target, nof_target]);
phase_mod_conj = exp(1i*2*pi*phase_mod/3);
phase_mod_conj = cat(5, ones(size(phase_mod_conj)), phase_mod_conj, conj(phase_mod_conj));

theta = reshape(theta, [1, nor_target, nof_target]);


rays = 1:nor_target;
kernel = zeros([nc*nsms, nc*nsms*kernel_length, kspace_segments, nor_target, nof_target], 'single');

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
kx_temp = kx(:, rays, iframe);
ky_temp = ky(:, rays, iframe);
N = NUFFT.init(repmat(kx_temp, [1, 1, nof_calib]), repmat(ky_temp, [1, 1, nof_calib]), 1, [6, 6], sx, sx);

kSpace_sb = NUFFT.NUFFT(im_sb, N);
phase_mod_temp = conj(phase_mod_conj(:, rays, iframe, :, :));
phase_idx_temp = phase_mod(:, rays, iframe);
kSpace_mb = sum(kSpace_sb .* phase_mod_temp, 5);
kSpace_mb_sl1 = kSpace_sb(:,:,:,:,1);
kSpace_mb_sl2 = kSpace_sb(:,:,:,:,2) .* phase_mod_temp(:,:,:,:,2);
kSpace_mb_sl3 = kSpace_sb(:,:,:,:,3) .* phase_mod_temp(:,:,:,:,3);

theta_temp = theta(:, rays, iframe);

for iray = 1:nor_target
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
    
    kSpace_source_phase_mod_12 = kSpace_mb_sl1(:, ray_idx_phase_mod_1, :, :);
    kSpace_source_phase_mod_22 = kSpace_mb_sl1(:, ray_idx_phase_mod_2, :, :);
    
    kSpace_source_phase_mod_13 = kSpace_mb_sl2(:, ray_idx_phase_mod_1, :, :);
    kSpace_source_phase_mod_23 = kSpace_mb_sl2(:, ray_idx_phase_mod_2, :, :);
    
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
    source_2 = cat(2, kSpace_mb_sl1(:, iray, :, :, :), kSpace_source_phase_mod_12, kSpace_source_phase_mod_22);
    source_3 = cat(2, kSpace_mb_sl2(:, iray, :, :, :), kSpace_source_phase_mod_13, kSpace_source_phase_mod_23);
    target = kSpace_sb(:, iray, :, :, :);
    
    
    for iseg = 1:kspace_segments
        segment_idx = (1:sx/kspace_segments) + (iseg-1)*sx/kspace_segments;
        if iseg ~= 1
            segment_idx = [segment_idx(1) - (1:point_extend), segment_idx];
        end
        if iseg ~= kspace_segments
            segment_idx = [segment_idx, segment_idx(end) + (1:point_extend)];
        end
        
        source_temp = source(segment_idx, :, :, :);
        source_temp_2 = source_2(segment_idx, :, :, :);
        source_temp_3 = source_3(segment_idx, :, :, :);
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
                source_kernel_temp_2 = source_temp_2(point_idx, :, itime, :);
                source_kernel_temp_3 = source_temp_3(point_idx, :, itime, :);
                source_kernel(:, (itime-1)*nkernel/nof_calib + ipatch) = source_kernel_temp(:);
%                 source_kernel(:, (itime-1)*nkernel/nof_calib + ipatch + nkernel) = source_kernel_temp_2(:)*1;
%                 source_kernel(:, (itime-1)*nkernel/nof_calib + ipatch + nkernel*2) = source_kernel_temp_3(:)*1;
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
kernel = reshape(kernel, [nc, nsms, nc*nsms*kernel_length, kspace_segments, nor_target, nof_target]);
kernel = permute(kernel, [1, 3, 4, 5, 2, 6]);
t1 = toc;
fprintf(repmat('\b',1,linelength));
fprintf(sprintf('calibration done, time = %.2f s\n', t1));
fprintf([repmat('-', [1, 75]), '\n'])


if 0

test_frame = 1:100;
ntest_frame = length(test_frame);
kSpace_temp = zeros([sx, nor_target, ntest_frame, nc, nsms]);
test_kspace = reshape(test_kspace, [sx, nor_target, nof_target, nc]);
kSpace_mb_test = test_kspace(:, :, test_frame, :);

for iframe = 1:ntest_frame
fprintf(sprintf('apply kernel to frame %g / %g \n', iframe, ntest_frame))
for islice = 1:nsms
for iray = 1:nor_target
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
    
    kSpace_source_phase_mod_1 = kSpace_mb_test(:, ray_idx_phase_mod_1, iframe, :);
    kSpace_source_phase_mod_2 = kSpace_mb_test(:, ray_idx_phase_mod_2, iframe, :);
    
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
    
    
    source = cat(2, kSpace_mb_test(:, iray, iframe, :), kSpace_source_phase_mod_1, kSpace_source_phase_mod_2);
%     target = kSpace_sb(:, i, :, :, islice);
    
    
    for iseg = 1:kspace_segments
        segment_idx = (1:sx/kspace_segments) + (iseg-1)*sx/kspace_segments;
        if iseg ~= 1
            segment_idx = [segment_idx(1) - point_extend, segment_idx];
        end
        if iseg ~= kspace_segments
            segment_idx = [segment_idx, segment_idx(end) + point_extend];
        end
        source_temp = source(segment_idx, :, :, :);
%         target_temp = target(segment_idx, :, :, :);
        nkernel = (length(segment_idx) - kernel_length + 1);
        
        % the nsms here is also number of calibration rays
        source_kernel = zeros([kernel_length * nc * nsms, nkernel], 'single');
%         target_kernel = zeros([nc, nkernel], 'single');
        
        for ipatch = 1:nkernel
            point_idx = (1:kernel_length) + (ipatch-1);
%             center = point_idx((kernel_length + 1)/2);
            source_kernel_temp = source_temp(point_idx, :, :, :);
            source_kernel(:, ipatch) = source_kernel_temp(:);
%             target_kernel(:, (itime-1)*nkernel/nof + ipatch) = target_temp(center, :, itime, :);
        end
        kspace_temp = kernel(:, :, iseg, iray, islice, test_frame(iframe)) * source_kernel;
        kSpace_temp(segment_idx(point_extend+1:end-point_extend), iray, iframe, :, islice) = kspace_temp.';
    end
end
end
end

N = NUFFT.init(kx(:, :, test_frame), ky(:, :, test_frame), 1, [6, 6], sx, sx);
im_test = NUFFT.NUFFT_adj(kSpace_temp, N);
% im_sb_test = NUFFT.NUFFT_adj(kSpace_sb(:,:,test_frame,:,:), N);
im_mb = NUFFT.NUFFT_adj(kSpace_mb_test .* phase_mod_conj(:, rays, test_frame, :, :), N);

% figure
% imagesc([sos(im_mb); sos(im_test); sos(im_sb_test)])
% axis image
% axis off
% colormap gray
% brighten(0.4)

im_test = sum(im_test .* conj(sens), 4);
% im_sb_test = sum(im_sb_test .* conj(sens), 4);
im_mb = sum(im_mb .* conj(sens), 4);

% figure
% imagesc(abs([im_mb(:,:); im_test(:,:); im_sb_test(:,:)]));
% axis image
% axis off
% colormap gray
% brighten(0.4)

figure
imagesc(abs([im_mb(:,:); im_test(:,:); im_sb_test(:,:)]));
axis image
axis off
colormap gray
brighten(0.4)
end
