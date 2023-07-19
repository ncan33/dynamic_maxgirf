function kSpace_cart = train_radial_ssg_traj_sp_cart_SING_op_dic(im_sb, theta, kx, ky, phase_mod_idx, kSpace_mb)

fprintf([repmat('-', [1, 75]), '\n'])
fprintf(sprintf('begin trajectory specified radial slice GRAPPA calibration: \n'));
fprintf([repmat('-', [1, 75]), '\n'])

kernel_length = 5;
nrays = 5;
max_d = 5;
lambda = 0.01;
kb_size = [2, 2];

% size
[sx, ~, nof_calib, nc, nsms] = size(im_sb);
kspace_segments = sx;
nor_target = size(kx, 2);
nof_target = size(kx, 3);

% single band cartesian k-space
N = NUFFT.init(kx, ky, 1, kb_size, sx, sx);
kSpace_cart_calib = fft2(im_sb.*N.Apodizer); clear im_sb

% calibration area
calib_begin_x = -10;
calib_begin_y = -10;

% shift kernel
shift_kx = 0:2:20;
shift_ky = 0:2:20;
[shift_kx, shift_ky] = meshgrid(shift_kx, shift_ky);
shift_kx = shift_kx(:)';
shift_ky = shift_ky(:)';
nshift = length(shift_kx);
% nshift = 73;

% phase modulation
phase_mod_idx = reshape(phase_mod_idx, [1, nor_target, nof_target]);
phase_mod_ = exp(-1i*2*pi*phase_mod_idx/3);
phase_mod_ = cat(5, ones(size(phase_mod_)), phase_mod_, conj(phase_mod_));

theta = reshape(theta, [1, nor_target, nof_target]);

% get which cartesian k-space to estimate
kx_cart = round(kx + sx/2 + 1);
ky_cart = round(ky + sx/2 + 1);
kx_cart(kx_cart < 1) = sx - kx_cart(kx_cart < 1);
ky_cart(ky_cart < 1) = sx - ky_cart(ky_cart < 1);
kx_cart(kx_cart > sx) = kx_cart(kx_cart > sx) - sx;
ky_cart(ky_cart > sx) = ky_cart(ky_cart > sx) - sx;

kx_cart = kx_cart - sx/2 - 1;
ky_cart = ky_cart - sx/2 - 1;

% init
kSpace_cart = zeros([sx, sx, nof_target, nc, nsms], 'like', kSpace_cart_calib);

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
    phase_idx_temp = phase_mod_idx(:, :, iframe);
    theta_temp = theta(:, :, iframe);
    
    kernel = cell(sx, nor_target);
    
    for iray = 1:nor_target
        fprintf(sprintf('ray = %g/%g ', iray, nor_target))
        tic
        
        % find n rays that are mostly close
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
        
        % init some cell
        n_point = zeros([sx, 1], 'single');
        kx_source_all = cell(sx, 1);
        ky_source_all = cell(sx, 1);
        phase_mod_source_all = cell(sx, 1);
        
        for ipoint = 1:sx
            kx_target = kx_cart_temp(ipoint, iray);
            ky_target = ky_cart_temp(ipoint, iray);
            
            d = sqrt((kx_target - kx_temp(:, ray_idx)).^2 + (ky_target - ky_temp(:, ray_idx)).^2);
            
            [~, point_idx] = mink(d, kernel_length);
            
            kx_source = zeros([kernel_length, nrays], 'like', kSpace_cart_calib);
            ky_source = zeros([kernel_length, nrays], 'like', kSpace_cart_calib);
            d_ = zeros([kernel_length, nrays], 'single');
            for ii = 1:nrays
                kx_source(:, ii) = kx_temp(point_idx(:, ii), ray_idx(ii));
                ky_source(:, ii) = ky_temp(point_idx(:, ii), ray_idx(ii));
                d_(:, ii) = d(point_idx(:, ii), ii);
            end
            
            idx_drop = d_ > max_d; extend = 1;
            while sum(idx_drop(:)) == kernel_length * nrays
                idx_drop = d_ > max_d + 0.1 * extend;
                extend = extend + 1;
            end
            n_point(ipoint) = sum(~idx_drop(:)) + 1;
            
            kx_source(idx_drop) = [];
            ky_source(idx_drop) = [];
            
            phase_mod_source = squeeze(repmat(phase_mod_temp(:, ray_idx, :, :, :), [kernel_length, 1]));
            phase_mod_source(repmat(idx_drop, [1, 1, nsms])) = [];
            phase_mod_source_all{ipoint} = reshape(phase_mod_source, [n_point(ipoint) - 1, 1, 1, 1, nsms]);
            
            kx_all = [kx_target; kx_source(:)];
            ky_all = [ky_target; ky_source(:)];
            
            % put the kernel to bottow left cornor in calibration area
            % and add all the shifts
%             kx_all = kx_all - min(kx_all) + calib_begin_x + shift_kx;
%             ky_all = ky_all - min(ky_all) + calib_begin_y + shift_ky;
            kx_all = kx_all + calib_begin_x + shift_kx;
            ky_all = ky_all + calib_begin_y + shift_ky;

%             r = sqrt(kx_target^2 + ky_target^2);
%             ringx = r * cos(0:0.05:2*pi) - kx_target;
%             ringy = r * sin(0:0.05:2*pi) - ky_target;
%             kx_all = kx_all + ringx;
%             ky_all = ky_all + ringy;
            
%             if ipoint < 72
%                 kx_all = kx_all + cos(theta_this) * ((1:73) - ipoint);
%                 ky_all = ky_all + sin(theta_this) * ((1:73) - ipoint);
%             elseif ipoint > 216
%                 kx_all = kx_all + cos(theta_this) * ((-72:0) + 288 - ipoint);
%                 ky_all = ky_all + sin(theta_this) * ((-72:0) + 288 - ipoint);
%             else
%                 kx_all = kx_all + cos(theta_this) * (-36:36);
%                 ky_all = ky_all + sin(theta_this) * (-36:36);
%             end
            
%             kx_all = kx_all + reshape([0,5,10], [1,1,3]);
%             ky_all = ky_all + reshape([0,5,10], [1,1,3]);
           
            kx_source_all{ipoint} = kx_all(:);
            ky_source_all{ipoint} = ky_all(:);
        end
        
        % interpolate calibration data
        n_point_all = sum(n_point) * nshift;
        kspace_calib = zeros([n_point_all, nc, nof_calib, nsms], 'like', kSpace_cart_calib);
        
        N = NUFFT.init(cat(1, kx_source_all{:}), cat(1, ky_source_all{:}), 1, kb_size, sx ,sx);
        for itime = 1:nof_calib
            kspace_calib(:, :, itime, :, :) = NUFFT.cart2rad(kSpace_cart_calib(:, :, itime, :, :), N);
        end
        
        % estimate GRAPPA kernel
%         tic
        for ipoint = 1:sx
            idx_temp = sum(n_point(1:ipoint-1)) * nshift + 1 : sum(n_point(1:ipoint)) * nshift;
            
            kspace_calib_temp = kspace_calib(idx_temp, :, :, :);
            kspace_calib_temp = reshape(kspace_calib_temp, [n_point(ipoint), nshift, nc, nof_calib, nsms]);
            
            source = kspace_calib_temp(2:end, :, :, :, :);
            source = sum(source .* phase_mod_source_all{ipoint}, 5);
            source = permute(source, [2, 4, 1, 3]);
            source = reshape(source, [nshift*nof_calib, (n_point(ipoint) - 1) * nc]);
            scale = max(abs(source(:)));
            source = source/scale;
            
            target = kspace_calib_temp(1, :, :, :, :);
            target = permute(target, [2, 4, 3, 5, 1]);
            target = reshape(target, [nshift*nof_calib, nc*nsms])/scale;
            
%             kernel{ipoint, iray} = target * pinv(source);
%             kernel{ipoint, iray} = lsqminnorm(source, target).';
            % Tikhonov regularizaiton
            kernel{ipoint, iray} = ((source' * source + lambda * eye((n_point(ipoint)-1)*nc))^-1 * source' * target).';
            
        end
        toc
    end
    [kSpace_cart(:, :, iframe, :, :)] = apply_radial_ssg_traj_specif_cart_SING_op(kSpace_mb(:, :, iframe, :, :), theta(:, :, iframe), kx(:, :, iframe), ky(:, :, iframe), kernel, max_d, kernel_length, nrays, phase_mod_idx);
    
end
t1 = toc;
fprintf(repmat('\b',1,linelength));
fprintf(sprintf('calibration done, time = %.2f s\n', t1));
fprintf([repmat('-', [1, 75]), '\n'])

% display test image
if 0
    para.Recon.kSpace_size = [sx, sx];
    para.Recon.nor = nor_target;
    para.over_sampling = 1;
    para.Recon.sx = sx;
    
    w = ramp_filter_for_pre_interp(para);
    
    test_im = fftshift2(ifft2(fftshift2(kSpace_cart).*w));
    test_im = squeeze(sos(test_im, 4));
    
    figure
    imagesc(test_im(:,:))
    colormap gray
    axis image
    axis off
    brighten(0.4)
end