function kernel = train_radial_ssg_trajectory_specified_to_cartesian_grid_SING_op(im_sb, theta, kx, ky, phase_mod_idx, k_test)

fprintf([repmat('-', [1, 75]), '\n'])
fprintf(sprintf('begin trajectory specified radial slice GRAPPA calibration: \n'));
fprintf([repmat('-', [1, 75]), '\n'])

tic

kernel_length = 5;
nrays = 5;
max_d = 144;

% calibration area
calib_begin_x = -30;
calib_begin_y = -30;

% shift kernel
shift_kx = 0:5:60;
shift_ky = 0:5:60;
[shift_kx, shift_ky] = meshgrid(shift_kx, shift_ky);
shift_kx = shift_kx(:)';
shift_ky = shift_ky(:)';
nshift = length(shift_kx);

[sx, ~, nof_calib, nc, nsms] = size(im_sb);
kspace_segments = sx;
nor_target = size(kx, 2);
nof_target = size(kx, 3);

% phase modulation
phase_mod_idx = reshape(phase_mod_idx, [1, nor_target, nof_target]);
phase_mod_ = exp(-1i*2*pi*phase_mod_idx/3);
phase_mod_ = cat(5, ones(size(phase_mod_)), phase_mod_, conj(phase_mod_));

theta = reshape(theta, [1, nor_target, nof_target]);

kernel = cell(sx, nor_target, nof_target);
% idx_source = cell(sx, nor_target, nof_target);

% get which cartesian k-space to estimate
kx_cart = round(kx + sx/2+1);
ky_cart = round(ky + sx/2+1);
kx_cart(kx_cart < 1) = sx - kx_cart(kx_cart < 1);
ky_cart(ky_cart < 1) = sx - ky_cart(ky_cart < 1);
kx_cart(kx_cart > sx) = kx_cart(kx_cart > sx) - sx;
ky_cart(ky_cart > sx) = ky_cart(ky_cart > sx) - sx;

kx_cart = kx_cart - sx/2 - 1;
ky_cart = ky_cart - sx/2 - 1;

% single band cartesian k-space
N = NUFFT.init(kx, ky, 1, [2, 2], sx, sx);
kSpace_cart = fft2(im_sb.*N.Apodizer); clear im_sb

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
    
    kx_temp = kx(:, :, iframe);
    ky_temp = ky(:, :, iframe);
    
    kx_cart_temp = kx_cart(:, :, iframe);
    ky_cart_temp = ky_cart(:, :, iframe);
    
    phase_mod_temp = phase_mod_(:, :, iframe, :, :);
    % phase_idx_temp = phase_mod(:, :, iframe);
    theta_temp = theta(:, :, iframe);
    
    for iray = 1:nor_target
        fprintf(sprintf('ray = %g/%g', iray, nor_target))
        tic
        
        % find 3 rays that are mostly close
        theta_this = theta_temp(iray);
        dtheta = mod(theta_this - theta_temp, pi);
        dtheta = min(dtheta, pi - dtheta);
        [~, ray_idx] = mink(dtheta, nrays);
        
        parfor ipoint = 1:sx
            kx_target = kx_cart_temp(ipoint, iray);
            ky_target = ky_cart_temp(ipoint, iray);
            
            d = sqrt((kx_target - kx_temp(:, ray_idx)).^2 + (ky_target - ky_temp(:, ray_idx)).^2);
            
            [~, point_idx] = mink(d, kernel_length);
            
            kx_source = zeros([kernel_length, nrays], 'like', kSpace_cart);
            ky_source = zeros([kernel_length, nrays], 'like', kSpace_cart);
            d_ = zeros([kernel_length, nrays], 'single');
            for ii = 1:nrays
                kx_source(:, ii) = kx_temp(point_idx(:, ii), ray_idx(ii));
                ky_source(:, ii) = ky_temp(point_idx(:, ii), ray_idx(ii));
                d_(:, ii) = d(point_idx(:, ii), ii);
            end
            
            idx_drop = d_ > max_d;
            nsource_point = sum(~idx_drop(:));
            
            kx_source(idx_drop) = [];
            ky_source(idx_drop) = [];
            
            phase_mod_source = repmat(phase_mod_temp(:, ray_idx, :, :, :), [kernel_length, 1]);
            phase_mod_source(repmat(idx_drop, [1, 1, nsms])) = [];
            phase_mod_source = reshape(phase_mod_source, [nsource_point, 1, 1, 1, nsms]);
            
            kx_all = [kx_target; kx_source(:)];
            ky_all = [ky_target; ky_source(:)];
            
            % put the kernel to bottow left cornor in calibration area
            kx_all = kx_all - min(kx_all) + calib_begin_x;
            ky_all = ky_all - min(ky_all) + calib_begin_y;
            
            kx_all = kx_all + shift_kx;
            ky_all = ky_all + shift_ky;

            tic
            N = NUFFT.init(kx_all(:), ky_all(:), 1, [2, 2], sx, sx);
            kspace_calib = zeros([(nsource_point+1)*nshift, nc, nof_calib, nsms], 'like', kSpace_cart);
            for itime = 1:nof_calib
                kspace_calib(:, :, itime, :, :) = NUFFT.cart2rad(kSpace_cart(:, :, itime, :, :), N);
            end
            toc
  
            kspace_calib = reshape(kspace_calib, [nsource_point+1, nshift, nc, nof_calib, nsms]);
            source = kspace_calib(2:end, :, :, :, :);
            target = kspace_calib(1, :, :, :, :);
            
            source = sum(source .* phase_mod_source, 5);
            
            source = permute(source, [1, 3, 2, 4]);
            source = reshape(source, [nsource_point*nc, nshift*nof_calib]);
            target = permute(target, [3, 5, 2, 4, 1]);
            target = reshape(target, [nc*nsms, nshift*nof_calib]);
            
            kernel{ipoint, iray, iframe} = target * pinv(source);
            % Tikhonov regularizaiton
            %         lambda = 1e13;
            %         kernel{ipoint, iray, iframe} = ((source.')' * source.' + lambda * eye(nsource_point*nc))^-1 * (source.')' * target.';
            %
            %         [point_idx, ray_idx_] = meshgrid(point_idx, ray_idx);
            %         idx_source{ipoint, iray, iframe} = sub2ind([sx, nor_target], point_idx, ray_idx_);
            %         idx_source{ipoint, iray, iframe}(idx_drop) = [];
        end
        toc
    end
    
end
t1 = toc;
fprintf(repmat('\b',1,linelength));
fprintf(sprintf('calibration done, time = %.2f s\n', t1));
fprintf([repmat('-', [1, 75]), '\n'])

% test calibration
temp_cart = apply_radial_ssg_traj_specif_cart_SING_op(k_test, theta, kx, ky, kernel, max_d, kernel_length, nrays);


para.Recon.kSpace_size = [sx, sx];
para.Recon.nor = nor_target;
para.over_sampling = 1;
para.Recon.sx = sx;

w = ramp_filter_for_pre_interp(para);

test_im = fftshift2(ifft2(fftshift2(temp_cart).*w));
test_im = squeeze(sos(test_im, 4));

figure
imagesc(test_im(:,:))
colormap gray
axis image
axis off
brighten(0.4)



keyboard