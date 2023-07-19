function kernel = train_radial_ssg_trajectory_specified_to_cartesian_grid_SING(im_sb, theta, kx, ky, phase_mod, k_test)

fprintf([repmat('-', [1, 75]), '\n'])
fprintf(sprintf('begin trajectory specified radial slice GRAPPA calibration: \n'));
fprintf([repmat('-', [1, 75]), '\n'])
                                                                           
tic

kernel_length = 3;

point_extend = (kernel_length - 1)/2;

[sx, ~, nof_calib, nc, nsms] = size(im_sb);
kspace_segments = sx;
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

kx_cart = kx_cart - sx/2 - 1;
ky_cart = ky_cart - sx/2 - 1;


% get Cartesian k-space
kSpace_cart_sb = fftshift2(fft2(fftshift2(im_sb)));

% display some text
text_display = sprintf('number of time frames: %16.0f\nnumber of calibration time frames: %4.0f\nnumber of coils: %22.0f\nnumber of slices: %21.0f\n\n', nof_target, nof_calib, nc, nsms);
fprintf(text_display)
text_display = sprintf('kernel length: %24.0f\nnumber of k-space segments: %11.0f\n\n', kernel_length, kspace_segments);
fprintf(text_display)

kx_shift = -144:10:144;
ky_shift = -144:3:144;
[kx_shift, ky_shift] = meshgrid(kx_shift, ky_shift);
kx_shift = kx_shift(:);
ky_shift = ky_shift(:);
nshift = length(kx_shift);
kx_shift = reshape(kx_shift, [1, nshift]);
ky_shift = reshape(ky_shift, [1, nshift]);

for iframe = 1:nof_target
if iframe ~= 1
    fprintf(repmat('\b',1,linelength));
end
linelength = fprintf(sprintf('calibrate frame %g / %g \n', iframe, nof_target));

% get trajectory specified calibration data
% interpolate single band image into target trajectory
kx_temp = kx(:, :, iframe);
ky_temp = ky(:, :, iframe);
% N = NUFFT.init(kx_temp, ky_temp, 1, [6, 6], sx, sx);
% 
% kSpace_sb = zeros([sx, nor_target, nof_calib, nc, nsms], 'like', im_sb);
% for itemp = 1:nof_calib
%     kSpace_sb(:, :, itemp, :, :) = NUFFT.NUFFT(im_sb(:, :, itemp, :, :), N);
% end

% cartesian target k-space data
% idx_cart = sub2ind([sx, sx], kx_cart(:, :, iframe), ky_cart(:, :, iframe));
% idx_cart = repmat(idx_cart, [1, 1, nof_calib, nc, nsms]);
% idx_cart = idx_cart + reshape((0:nof_calib-1)*sx*sx, [1, 1, nof_calib]);
% idx_cart = idx_cart + reshape((0:nc-1)*sx*sx*nof_calib, [1, 1, 1, nc]);
% idx_cart = idx_cart + reshape((0:nsms-1)*sx*sx*nof_calib*nc, [1, 1, 1, 1, nsms]);
% kSpace_cart_target = kSpace_cart_sb(idx_cart);

kx_cart_temp = kx_cart(:, :, iframe);
ky_cart_temp = ky_cart(:, :, iframe);

phase_mod_temp = conj(phase_mod_conj(:, :, iframe, :, :));
phase_idx_temp = phase_mod(:, :, iframe);
% kSpace_mb = sum(kSpace_sb .* phase_mod_temp, 5);
theta_temp = theta(:, :, iframe);

for iray = 1:nor_target
    iray
    tic
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
    
    for ipoint = 1:sx

        kx_target = kx_cart_temp(ipoint, iray);
        ky_target = ky_cart_temp(ipoint, iray);
        
        kx_source_phase_mod_1 = kx_temp(:, ray_idx_phase_mod_1);
        ky_source_phase_mod_1 = ky_temp(:, ray_idx_phase_mod_1);
        
        kx_source_phase_mod_2 = kx_temp(:, ray_idx_phase_mod_2);
        ky_source_phase_mod_2 = ky_temp(:, ray_idx_phase_mod_2);
        
        dx_phase_mod_1 = kx_target - kx_source_phase_mod_1;
        dy_phase_mod_1 = ky_target - ky_source_phase_mod_1;
        
        d_phase_mod_1 = sqrt(dx_phase_mod_1.^2 + dy_phase_mod_1.^2);
        
        [~, point_idx_phase_mod_1] = mink(d_phase_mod_1, kernel_length);
        point_idx_phase_mod_1 = sort(point_idx_phase_mod_1);
        
        dx_phase_mod_2 = kx_target - kx_source_phase_mod_2;
        dy_phase_mod_2 = ky_target - ky_source_phase_mod_2;
        
        d_phase_mod_2 = sqrt(dx_phase_mod_2.^2 + dy_phase_mod_2.^2);
        
        [~, point_idx_phase_mod_2] = mink(d_phase_mod_2, kernel_length);
        point_idx_phase_mod_2 = sort(point_idx_phase_mod_2);
        
        kx_phase_mod_1 = kx_source_phase_mod_1(point_idx_phase_mod_1);
        ky_phase_mod_1 = ky_source_phase_mod_1(point_idx_phase_mod_1);
        
        kx_phase_mod_2 = kx_source_phase_mod_2(point_idx_phase_mod_2);
        ky_phase_mod_2 = ky_source_phase_mod_2(point_idx_phase_mod_2);
        
        if ipoint <= point_extend
            point_idx_phase_mod_0 = 1:kernel_length;
        elseif ipoint > sx - point_extend
            point_idx_phase_mod_0 = sx - kernel_length + 1 : sx;
        else
            point_idx_phase_mod_0 = ipoint - point_extend : ipoint + point_extend;
        end
        kx_phase_mod_0 = kx_temp(point_idx_phase_mod_0, iray);
        ky_phase_mod_0 = ky_temp(point_idx_phase_mod_0, iray);
        
        kx_source = cat(2, kx_phase_mod_0, kx_phase_mod_1, kx_phase_mod_2);
        ky_source = cat(2, ky_phase_mod_0, ky_phase_mod_1, ky_phase_mod_2);
        
        kx_all = [kx_target; kx_source(:)];
        ky_all = [ky_target; ky_source(:)];
        
        kx_all = kx_all + kx_shift;
        ky_all = ky_all + ky_shift;
        
        N = NUFFT.init(kx_all(:), ky_all(:), 1, [6, 6], sx, sx);
        kspace_calib = zeros([(kernel_length*nsms+1)*nshift, nc, nof_calib, nsms], 'like', im_sb);
        for itime = 1:nof_calib
            kspace_calib(:, :, itime, :, :) = NUFFT.NUFFT(im_sb(:, :, itime, :, :), N);
        end
        
        kspace_calib = reshape(kspace_calib, [kernel_length*nsms+1, nshift, nc, nof_calib, nsms]);
        source = kspace_calib(2:end, :, :, :, :);
        target = kspace_calib(1, :, :, :, :);
        
        phase_mod_temp_source = phase_mod_temp(:, [iray, ray_idx_phase_mod_1, ray_idx_phase_mod_2], :, :, :);
        phase_mod_temp_source = repmat(phase_mod_temp_source, [kernel_length, 1, 1, 1, 1]);
        phase_mod_temp_source = reshape(phase_mod_temp_source, [kernel_length*nsms, 1, 1, 1, nsms]);
        
        source = sum(source .* phase_mod_temp_source, 5);
        source = permute(source, [1, 3, 2, 4]);
        source = reshape(source, [kernel_length*nsms*nc, nshift*nof_calib]);
        target = permute(target, [3, 5, 2, 4, 1]);
        target = reshape(target, [nc*nsms, nshift*nof_calib]);
        
        kernel(:, :, ipoint, iray) = target * pinv(source);
        
    end
    toc
end

end
t1 = toc;
fprintf(repmat('\b',1,linelength));
fprintf(sprintf('calibration done, time = %.2f s\n', t1));
fprintf([repmat('-', [1, 75]), '\n'])



temp_cart = zeros([sx, sx, nof_target, nc, nsms], 'single');
temp_cart_1 = zeros([sx, sx], 'single');
k_temp = apply_radial_ssg_trajectroy_speficied_to_cartesian_grid(k_test, phase_mod, theta, kernel);
for iframe = 1:nof_target
    for ic = 1:nc
        for isms = 1:nsms
            temp_cart_1(idx_cart) = k_temp(:, :, 1, ic, isms);
            temp_cart(:, :, iframe, ic, isms) = temp_cart_1;
        end
    end
end

para.Recon.kSpace_size = [sx, sx];
para.Recon.nor = nor_target;
para.over_sampling = 1;
para.Recon.sx = sx;

w = ramp_filter_for_pre_interp(para);

test_im = fftshift2(ifft2(fftshift2(temp_cart).*w));



keyboard