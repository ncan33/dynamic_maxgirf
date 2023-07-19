function [temp_cart, kSpace_sb] = apply_radial_ssg_traj_specif_cart_SING_op(kSpace_mb, theta, kx, ky, kernel, max_d, kernel_length, nrays, phase_mod_idx)

[sx, nor, nof, nc] = size(kSpace_mb);
nsms = size(kernel{1}, 1)/nc;

kSpace_sb = zeros([sx, nor, nof, nc, nsms], 'like', kSpace_mb);

% max_d = 144;
% kernel_length = 5;
% nrays = 3;

kx_cart = round(kx + sx/2+1);
ky_cart = round(ky + sx/2+1);
kx_cart(kx_cart < 1) = sx - kx_cart(kx_cart < 1);
ky_cart(ky_cart < 1) = sx - ky_cart(ky_cart < 1);
kx_cart(kx_cart > sx) = kx_cart(kx_cart > sx) - sx;
ky_cart(ky_cart > sx) = ky_cart(ky_cart > sx) - sx;

kx_cart = kx_cart - sx/2 - 1;
ky_cart = ky_cart - sx/2 - 1;

for iframe = 1:nof
    
    % get trajectory specified calibration data
    % interpolate single band image into target trajectory
    kx_temp = kx(:, :, iframe);
    ky_temp = ky(:, :, iframe);

    theta_temp = theta(:, :, iframe);
    
    kx_cart_temp = kx_cart(:, :, iframe);
    ky_cart_temp = ky_cart(:, :, iframe);
    
    phase_idx_temp = phase_mod_idx(:, :, iframe);
    for iray = 1:nor
        
        % find 3 rays that are mostly close
        theta_this = theta_temp(iray);
        dtheta = mod(theta_this - theta_temp, pi);
        dtheta = min(dtheta, pi - dtheta);
        [~, ray_idx] = mink(dtheta, nrays);
        
        % find rays that are close and not the same phase mod
%         phase_idx_this = phase_idx_temp(iray);
%         ray_idx_1 = iray;
%         [~, ray_idx_2] = min(dtheta(phase_idx_temp ~= phase_idx_this));
%         temp = find(phase_idx_temp ~= phase_idx_this); ray_idx_2 = temp(ray_idx_2);
%         phase_idx_that = phase_idx_temp(ray_idx_2);
%         [~, ray_idx_3] = min(dtheta(phase_idx_temp ~= phase_idx_this & phase_idx_temp ~= phase_idx_that));
%         temp = find(phase_idx_temp ~= phase_idx_this & phase_idx_temp ~= phase_idx_that); ray_idx_3 = temp(ray_idx_3);
%         ray_idx = [ray_idx_1, ray_idx_2, ray_idx_3];
        
        
        for ipoint = 1:sx
            
            kx_target = kx_cart_temp(ipoint, iray);
            ky_target = ky_cart_temp(ipoint, iray);
%             kx_target = kx_temp(ipoint, iray);
%             ky_target = ky_temp(ipoint, iray);
            
            d = sqrt((kx_target - kx_temp(:, ray_idx)).^2 + (ky_target - ky_temp(:, ray_idx)).^2);
            
            [~, point_idx] = mink(d, kernel_length);
            
            kx_source = zeros([kernel_length, nrays], 'like', kSpace_mb);
            ky_source = zeros([kernel_length, nrays], 'like', kSpace_mb);
            d_ = zeros([kernel_length, nrays], 'single');
            source = zeros([kernel_length, nrays, nc], 'like', kSpace_mb);
            for ii = 1:nrays
                kx_source(:, ii) = kx_temp(point_idx(:, ii), ray_idx(ii));
                ky_source(:, ii) = ky_temp(point_idx(:, ii), ray_idx(ii));
                d_(:, ii) = d(point_idx(:, ii), ii);
                source(:, ii, :) = kSpace_mb(point_idx(:, ii), ray_idx(ii), iframe, :);
            end
            
            idx_drop = d_ > max_d; extend = 1;
            while sum(idx_drop(:)) == kernel_length * nrays
                idx_drop = d_ > max_d + 0.1 * extend;
                extend = extend + 1;
            end

            source(repmat(idx_drop, [1, 1, nc])) = [];

            temp = kernel{ipoint, iray, iframe} * source(:);
            
            kSpace_sb(ipoint, iray, iframe, :, :) = reshape(temp, [nc, nsms]);
        end
    end
    
end

temp_cart = zeros([sx, sx, nof, nc, nsms], 'single');
mask = zeros([sx, sx, nof, nc, nsms], 'single');

for iframe = 1:nof
    idx_cart = sub2ind([sx, sx], kx_cart(:, :, iframe) + sx/2+1, ky_cart(:, :, iframe) + sx/2+1);
    for ic = 1:nc
        for isms = 1:nsms
            for iray = 1:nor
                temp_cart_1 = zeros([sx, sx], 'single');
                mask_1 = temp_cart_1;
                temp_cart_1(idx_cart(:, iray)) = kSpace_sb(:, iray, iframe, ic, isms);
                mask_1(idx_cart(:, iray)) = 1;
                temp_cart(:, :, iframe, ic, isms) = temp_cart(:, :, iframe, ic, isms) + temp_cart_1;
                mask(:, :, iframe, ic, isms) = mask(:, :, iframe, ic, isms) + mask_1;
            end
        end
    end
end

mask(mask == 0) = 1;
temp_cart = temp_cart./mask;