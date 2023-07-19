function kernel = train_radial_ssg(kSpace_mb, kSpace_sb, theta, phase_mod, sens)

kernel_length = 5;
angle_range = 16;

ray_range = 1:288;
% ray_range = 73+32:216-32;

[sx, nor, nof, nc, nsms] = size(kSpace_mb);

kSpace_mb = sum(kSpace_mb, 5);
kSpace_mb = reshape(kSpace_mb, [sx, nor*nof, nc]);
kSpace_sb = reshape(kSpace_sb, [sx, nor*nof, nc, nsms]);
theta = reshape(theta, [1, nor*nof]);
phase_mod = reshape(phase_mod, [1, nor*nof]);

phase_mod_conj = exp(-1i*2*pi*phase_mod/3);
phase_mod_conj = cat(5, ones(size(phase_mod_conj)), phase_mod_conj, conj(phase_mod_conj));

kernel = zeros([nc, nc*kernel_length, nsms, angle_range, nsms], 'single');

for isms = 1:nsms
    idx_sms = phase_mod == isms-1;
    
    kSpace_mb_temp = kSpace_mb(ray_range, idx_sms, :);
    kSpace_sb_temp = kSpace_sb(ray_range, idx_sms, :, :);
    phase_mod_temp = phase_mod_conj(:, idx_sms, :, :, :);
    theta_temp = theta(idx_sms);
    
    for itheta = 1:angle_range
        theta_range = [0, pi/angle_range] + pi/angle_range * (itheta-1);
        idx_theta = theta_temp >= theta_range(1) & theta_temp < theta_range(2);
        
        kSpace_mb_temp_2 = kSpace_mb_temp(:, idx_theta, :);
        kSpace_sb_temp_2 = kSpace_sb_temp(:, idx_theta, :, :);
        phase_mod_temp_2 = phase_mod_temp(:, idx_theta, :, :, :);
%         theta_temp_2 = theta_temp(idx_theta);
        
        npatch = length(ray_range) - kernel_length + 1;
        
        nor_theta = sum(idx_theta);
        
        for islice = 1:nsms
            target = zeros([1, nor_theta, npatch, nc], 'single');
            source = zeros([kernel_length, nor_theta, npatch, nc], 'single');
            for ipatch = 1:npatch
                idx_patch = (1:kernel_length) + (ipatch-1);
                center = (kernel_length + 1)/2;
                source(:, :, ipatch, :) = kSpace_mb_temp_2(idx_patch, :, :);% .* phase_mod_temp_2(:, :, :, :, islice);
                target(1, :, ipatch, :) = kSpace_sb_temp_2(idx_patch(center), :, :, islice);
            end
            source = permute(source, [1, 4, 2, 3]);
            source = reshape(source, [kernel_length*nc, nor_theta*npatch]);
            target = permute(target, [1, 4, 2, 3]);
            target = reshape(target, [nc, nor_theta*npatch]);
            kernel(:, :, islice, itheta, isms) = target * pinv(source);
        end
    end
end
        



%% kernel apply

clear kSpace_sb_temp kSpace_sb_temp_2
for isms = 1:nsms
    idx_sms = phase_mod == isms-1;
    
    kSpace_mb_temp = kSpace_mb(:, idx_sms, :) .* phase_mod_conj(:,idx_sms,:,:,isms);
    phase_mod_temp = phase_mod_conj(:, idx_sms, :, :, :);
    theta_temp = theta(idx_sms);
    
    for itheta = 1:angle_range
        theta_range = [0, pi/angle_range] + pi/angle_range * (itheta-1);
        idx_theta = theta_temp >= theta_range(1) & theta_temp < theta_range(2);
        
        kSpace_mb_temp_2 = kSpace_mb_temp(:, idx_theta, :);
        phase_mod_temp_2 = phase_mod_temp(:, idx_theta, :, :, :);
        
        npatch = sx - kernel_length + 1;
        
        nor_theta = sum(idx_theta);
        
        for islice = 1:nsms
            source = zeros([kernel_length, nor_theta, npatch, nc], 'single');
            for ipatch = 1:npatch
                idx_patch = (1:kernel_length) + (ipatch-1);
                source(:, :, ipatch, :) = kSpace_mb_temp_2(idx_patch, :, :);% .* phase_mod_temp_2(:, :, :, :, islice);
            end
            kernel_temp = kernel(:, :, islice, itheta, isms);
            
            source = permute(source, [1, 4, 2, 3]);
            source = reshape(source, [kernel_length*nc, nor_theta*npatch]);
            
            target = kernel_temp * source;
            target = reshape(target, [nc, nor_theta, npatch]);
            target = permute(target, [2, 3, 1]);
            for ipatch = 1:npatch
                idx_patch = (1:kernel_length) + (ipatch-1);
                center = (kernel_length + 1)/2;                
                kSpace_sb_temp_2(idx_patch(center), :, :, islice) = target(:, ipatch, :);
            end
        end
        
        kSpace_sb_temp(:, idx_theta, :, :) = kSpace_sb_temp_2; clear kSpace_sb_temp_2
    end
    kSpace_sb_test(:, idx_sms, :, :) = kSpace_sb_temp; clear kSpace_sb_temp
end

kSpace_sb_test(end+1:end+1+(kernel_length-1)/2-1, :, :, :) = 0;

rays = 1201:1230;
[kx, ky] = get_k_coor(sx, theta(:,rays), 0, sx/2+1);
N = NUFFT.init(kx, ky, 1, [6, 6], sx, sx);

im_sb = NUFFT.NUFFT_adj(permute(kSpace_sb(:,rays,:,:), [1, 2, 5, 3, 4]), N);
im_mb_test = NUFFT.NUFFT_adj(permute(kSpace_mb(:,rays,:), [1, 2, 4, 3]).*phase_mod_conj(:, rays, :, :, :), N);

phase_mod_test = repmat(exp(1i*[0, 2*pi/3, 2*pi/3]), [1, 1000]);
im_sb_test = NUFFT.NUFFT_adj(permute(kSpace_sb_test(:, rays, :, :), [1, 2, 5, 3, 4]) .* phase_mod_test(:, rays, :, :, :), N);

im_mb_test = sum(im_mb_test .* conj(sens),4);
im_sb_test = sum(im_sb_test .* conj(sens),4);
im_sb = sum(im_sb .* conj(sens),4);

% im_mb_test = sos(im_mb_test, 4);
% im_sb_test = sos(im_sb_test, 4);
% im_sb = sos(im_sb, 4);

figure
imagesc(abs([im_mb_test(:,:); im_sb_test(:,:); im_mb_test(:,:) - im_sb_test(:,:)]))
axis image
colormap gray
brighten(0.4)
axis off
title(sprintf('%g ray/frame, calibration size = %g,\n k-space segments = %g, kernel size = %g', length(rays), length(ray_range), angle_range, kernel_length) , 'FontSize', 30)
set(gcf, 'Position', [0,0,900,900]);

