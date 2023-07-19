function [Data] = get_ssg_non_cart(Data, para)

[sx, nor_total, nc] = size(Data.kSpace);

if isfield(para.Recon.ssg, 'calib_size')
    calib_size = para.Recon.ssg.calib_size;
else
    calib_size = floor(sx/sqrt(2)/2)*2;
    calib_size = [calib_size, calib_size];
    para.Recon.ssg.calib_size = calib_size;
end

if ~isfield(para.Recon.ssg, 'patch_size')
    para.Recon.ssg.patch_size = [5, 5];
end

if ~isfield(para.Recon.ssg, 'alpha')
    para.Recon.ssg.patch_size = 4;
end

if ~isfield(para.Recon.ssg, 'section')
    para.Recon.ssg.section = 4;
end

if ~isfield(para.Recon.ssg, 'nor_calib')
    para.Recon.ssg.nor_calib = 4;
end

nsms = para.Recon.nSMS;

nor_calib = para.Recon.ssg.nor_calib;
nor_calib_begin = 1;
nor_calib_idx = nor_calib_begin:nor_calib_begin + nor_calib - 1;

nof_total = floor(nor_total/nor_calib);
nor_total = nof_total * nor_calib;

Data.kSpace(:, nor_total+1:end, :, :, :) = [];
Data.phase_mod(:, nor_total+1:end, :, :, :) = [];
Data.kx(:, nor_total+1:end) = [];
Data.ky(:, nor_total+1:end) = [];

Data.kSpace = reshape(Data.kSpace, [sx, nor_calib, nof_total, nc]);
Data.phase_mod = reshape(Data.phase_mod, [1, nor_calib, nof_total, 1, nsms]);
Data.kx = reshape(Data.kx, [sx, nor_calib, nof_total]);
Data.ky = reshape(Data.ky, [sx, nor_calib, nof_total]);

Data.N = NUFFT.init(Data.kx, Data.ky, 1, [6, 6], sx, sx);
Data.first_est = NUFFT.NUFFT_adj(Data.kSpace .* conj(Data.phase_mod), Data.N)/nsms;

%% set reconstruction parameters
% para.setting.ifplot = 1;
% para.setting.ifGPU = 0;
para.Recon.epsilon = eps('single');
para.Recon.step_size = 2;
para.Recon.noi = 30;
para.Recon.type = 'NUFFT coil';
para.Recon.break = 1;
para.Recon.weight_tTV = 0;
para.Recon.weight_sTV = 0;

%% reconstruct aliasing-free coil images
img = STCR_conjugate_gradient(Data, para);

%% get calibration k-space data
% single band radial k-space data
kSpace_raid_sb = NUFFT.NUFFT(img, Data.N);
% phase modulated single band radial k-space
kSpace_raid_pm = kSpace_raid_sb .* Data.phase_mod;

img_calib = zeros([sx, sx, nof_total, nc, nsms, nsms], 'single');

for i=1:nsms
    kSpace_raid_temp = kSpace_raid_pm .* conj(Data.phase_mod(:,:,:,i));
    img_calib(:,:,:,:,:,i) = NUFFT.NUFFT_adj(kSpace_raid_temp, Data.N);
end

% calibration k-space
kSpace_SSG_calib = fftshift2(fft2(fftshift2(img_calib)));
kSpace_calib_source = crop_half_FOV(kSpace_SSG_calib, calib_size);

para.Recon.ssg.imsize = [sx, sx];

switch para.Recon.ssg.section
    case 1
        [~, Data.ssg] = ssg_throught_time_train(kSpace_calib_source, para.Recon.ssg);
    case 4
        [~, Data.ssg] = ssg_throught_time_train_section(kSpace_calib_source, para.Recon.ssg);
    case 8
        [~, Data.ssg] = ssg_throught_time_train_8_section_2(kSpace_calib_source, para.Recon.ssg);
end


%% test kernel
if para.setting.ifplot && para.Recon.ssg.section == 4
    ssg = para.Recon.ssg;
    
    % nor = 300, calibration size [140, 140], alpha = 4
    nor = 30;
    N_test = NUFFT.init(Data.kx(:,1:nor,:), Data.ky(:,1:nor,:), 1, [6,6], sx, sx);
    img_test = NUFFT.NUFFT_adj(Data.kSpace(:,1:nor,:,:).*conj(Data.phase_mod(:,1:nor,:,:,:)), N_test);
    img_test_sens = squeeze(sum(img_test .* conj(Data.sens_map), 4));
    img_test = section_kernel_apply(img_test, Data.ssg);
    img_test_ssg_sens = squeeze(sum(img_test.*conj(Data.sens_map), 4));
    
    f = figure;
    imshow = abs([(img_test_sens(:,:,3,:)); (img_test_ssg_sens(:,:,3,:)); (abs(img_test_sens(:,:,3,:)- img_test_ssg_sens(:,:,3,:)))]);
    imshow = imshow(:,:);
    imagesc(imshow);axis image
    colormap gray; brighten(0.4); axis off
    title(sprintf('%g ray/frame, calibration size = [%g,%g],\n k-space segments = %g, \\alpha = %g, kernel size = [%g,%g]', nor, calib_size(1), calib_size(2), ssg.section, ssg.alpha, ssg.patch_size(1), ssg.patch_size(2)), 'FontSize', 30)
    set(gcf, 'Position', [0,0,900,900]);
%     hgexport(f,sprintf('./Figure/%g_ray_calibration_size_%g_%g_segments_%g_a_%g_kernel_size_%g_%g.eps', nor, calib_size(1), calib_size(2), ssg.section, ssg.alpha, ssg.patch_size(1), ssg.patch_size(2)))
    
    % nor = 300, calibration size [140, 140], alpha = 1
    ssg.alpha = 1;
    [~, Data.ssg] = ssg_throught_time_train_section(kSpace_calib_source, ssg);
    img_test = NUFFT.NUFFT_adj(Data.kSpace(:,1:nor,:,:).*conj(Data.phase_mod(:,1:nor,:,:,:)), N_test);
    img_test = section_kernel_apply(img_test, Data.ssg);
    img_test_ssg_sens = squeeze(sum(img_test.*conj(Data.sens_map), 4));
    
    f = figure;
    imshow = abs([fliplr(img_test_sens(:,:,1,:)); fliplr(img_test_ssg_sens(:,:,1,:)); fliplr(abs(img_test_sens(:,:,1,:)- img_test_ssg_sens(:,:,1,:)))]);
    imshow = imshow(:,:);
    imagesc(imshow);axis image
    colormap gray; brighten(0.4); axis off
    title(sprintf('%g ray/frame, calibration size = [%g,%g],\n k-space segments = %g, \\alpha = %g, kernel size = [%g,%g]', nor, calib_size(1), calib_size(2), ssg.section, ssg.alpha, ssg.patch_size(1), ssg.patch_size(2)), 'FontSize', 30)
    set(gcf, 'Position', [0,0,900,900]);
%     hgexport(f,sprintf('./Figure/%g_ray_calibration_size_%g_%g_segments_%g_a_%g_kernel_size_%g_%g.eps', nor, calib_size(1), calib_size(2), ssg.section, ssg.alpha, ssg.patch_size(1), ssg.patch_size(2)))

    % nor = 300, calibration size [140, 140], alpha = 10
    ssg.alpha = 10;
    [~, Data.ssg] = ssg_throught_time_train_section(kSpace_calib_source, ssg);
    img_test = NUFFT.NUFFT_adj(Data.kSpace(:,1:nor,:,:).*conj(Data.phase_mod(:,1:nor,:,:,:)), N_test);
    img_test = section_kernel_apply(img_test, Data.ssg);
    img_test_ssg_sens = squeeze(sum(img_test.*conj(Data.sens_map), 4));
    
    f = figure;
    imshow = abs([fliplr(img_test_sens(:,:,1,:)); fliplr(img_test_ssg_sens(:,:,1,:)); fliplr(abs(img_test_sens(:,:,1,:)- img_test_ssg_sens(:,:,1,:)))]);
    imshow = imshow(:,:);
    imagesc(imshow);axis image
    colormap gray; brighten(0.4); axis off
    title(sprintf('%g ray/frame, calibration size = [%g,%g],\n k-space segments = %g, \\alpha = %g, kernel size = [%g,%g]', nor, calib_size(1), calib_size(2), ssg.section, ssg.alpha, ssg.patch_size(1), ssg.patch_size(2)), 'FontSize', 30)
    set(gcf, 'Position', [0,0,900,900]);
%     hgexport(f,sprintf('./Figure/%g_ray_calibration_size_%g_%g_segments_%g_a_%g_kernel_size_%g_%g.eps', nor, calib_size(1), calib_size(2), ssg.section, ssg.alpha, ssg.patch_size(1), ssg.patch_size(2)))

    
    % nor = 300, calibration size [80, 80], alpha = 4
    calib_size = [80,80];
    kSpace_calib_source = crop_half_FOV(kSpace_SSG_calib, calib_size);
    ssg.alpha = 4;
    [~, Data.ssg] = ssg_throught_time_train_section(kSpace_calib_source, ssg);
    img_test = NUFFT.NUFFT_adj(Data.kSpace(:,1:nor,:,:).*conj(Data.phase_mod(:,1:nor,:,:,:)), N_test);
    img_test = section_kernel_apply(img_test, Data.ssg);
    img_test_ssg_sens = squeeze(sum(img_test.*conj(Data.sens_map), 4));
    
    f = figure;
    imshow = abs([fliplr(img_test_sens(:,:,1,:)); fliplr(img_test_ssg_sens(:,:,1,:)); fliplr(abs(img_test_sens(:,:,1,:)- img_test_ssg_sens(:,:,1,:)))]);
    imshow = imshow(:,:);
    imagesc(imshow);axis image
    colormap gray; brighten(0.4); axis off
    title(sprintf('%g ray/frame, calibration size = [%g,%g],\n k-space segments = %g, \\alpha = %g, kernel size = [%g,%g]', nor, calib_size(1), calib_size(2), ssg.section, ssg.alpha, ssg.patch_size(1), ssg.patch_size(2)), 'FontSize', 30)
    set(gcf, 'Position', [0,0,900,900]);
%     hgexport(f,sprintf('./Figure/%g_ray_calibration_size_%g_%g_segments_%g_a_%g_kernel_size_%g_%g.eps', nor, calib_size(1), calib_size(2), ssg.section, ssg.alpha, ssg.patch_size(1), ssg.patch_size(2)))

    % nor = 300, calibration size [80, 80], alpha = 4
    calib_size = [40,40];
    kSpace_calib_source = crop_half_FOV(kSpace_SSG_calib, calib_size);
    ssg.alpha = 4;
    [~, Data.ssg] = ssg_throught_time_train_section(kSpace_calib_source, ssg);
    img_test = NUFFT.NUFFT_adj(Data.kSpace(:,1:nor,:,:).*conj(Data.phase_mod(:,1:nor,:,:,:)), N_test);
    img_test = section_kernel_apply(img_test, Data.ssg);
    img_test_ssg_sens = squeeze(sum(img_test.*conj(Data.sens_map), 4));
    
    f = figure;
    imshow = abs([fliplr(img_test_sens(:,:,1,:)); fliplr(img_test_ssg_sens(:,:,1,:)); fliplr(abs(img_test_sens(:,:,1,:)- img_test_ssg_sens(:,:,1,:)))]);
    imshow = imshow(:,:);
    imagesc(imshow);axis image
    colormap gray; brighten(0.4); axis off
    title(sprintf('%g ray/frame, calibration size = [%g,%g],\n k-space segments = %g, \\alpha = %g, kernel size = [%g,%g]', nor, calib_size(1), calib_size(2), ssg.section, ssg.alpha, ssg.patch_size(1), ssg.patch_size(2)), 'FontSize', 30)
    set(gcf, 'Position', [0,0,900,900]);
%     hgexport(f,sprintf('./Figure/%g_ray_calibration_size_%g_%g_segments_%g_a_%g_kernel_size_%g_%g.eps', nor, calib_size(1), calib_size(2), ssg.section, ssg.alpha, ssg.patch_size(1), ssg.patch_size(2)))

    
    % nor = 300, calibration size [80, 80], alpha = 4
    ssg.patch_size = [3,3];
    calib_size = [140,140];
    kSpace_calib_source = crop_half_FOV(kSpace_SSG_calib, calib_size);
    ssg.alpha = 4;
    [~, Data.ssg] = ssg_throught_time_train_section(kSpace_calib_source, ssg);
    img_test = NUFFT.NUFFT_adj(Data.kSpace(:,1:nor,:,:).*conj(Data.phase_mod(:,1:nor,:,:,:)), N_test);
    img_test = section_kernel_apply(img_test, Data.ssg);
    img_test_ssg_sens = squeeze(sum(img_test.*conj(Data.sens_map), 4));
    
    f = figure;
    imshow = abs([fliplr(img_test_sens(:,:,1,:)); fliplr(img_test_ssg_sens(:,:,1,:)); fliplr(abs(img_test_sens(:,:,1,:)- img_test_ssg_sens(:,:,1,:)))]);
    imshow = imshow(:,:);
    imagesc(imshow);axis image
    colormap gray; brighten(0.4); axis off
    title(sprintf('%g ray/frame, calibration size = [%g,%g],\n k-space segments = %g, \\alpha = %g, kernel size = [%g,%g]', nor, calib_size(1), calib_size(2), ssg.section, ssg.alpha, ssg.patch_size(1), ssg.patch_size(2)), 'FontSize', 30)
    set(gcf, 'Position', [0,0,900,900]);
%     hgexport(f,sprintf('./Figure/%g_ray_calibration_size_%g_%g_segments_%g_a_%g_kernel_size_%g_%g.eps', nor, calib_size(1), calib_size(2), ssg.section, ssg.alpha, ssg.patch_size(1), ssg.patch_size(2)))

    
    % nor = 300, calibration size [80, 80], alpha = 4
    ssg.patch_size = [9,9];
    calib_size = [140,140];
    kSpace_calib_source = crop_half_FOV(kSpace_SSG_calib, calib_size);
    ssg.alpha = 4;
    [~, Data.ssg] = ssg_throught_time_train_section(kSpace_calib_source, ssg);
    img_test = NUFFT.NUFFT_adj(Data.kSpace(:,1:nor,:,:).*conj(Data.phase_mod(:,1:nor,:,:,:)), N_test);
    img_test = section_kernel_apply(img_test, Data.ssg);
    img_test_ssg_sens = squeeze(sum(img_test.*conj(Data.sens_map), 4));
    
    f = figure;
    imshow = abs([fliplr(img_test_sens(:,:,1,:)); fliplr(img_test_ssg_sens(:,:,1,:)); fliplr(abs(img_test_sens(:,:,1,:)- img_test_ssg_sens(:,:,1,:)))]);
    imshow = imshow(:,:);
    imagesc(imshow);axis image
    colormap gray; brighten(0.4); axis off
    title(sprintf('%g ray/frame, calibration size = [%g,%g],\n k-space segments = %g, \\alpha = %g, kernel size = [%g,%g]', nor, calib_size(1), calib_size(2), ssg.section, ssg.alpha, ssg.patch_size(1), ssg.patch_size(2)), 'FontSize', 30)
    set(gcf, 'Position', [0,0,900,900]);
%     hgexport(f,sprintf('./Figure/%g_ray_calibration_size_%g_%g_segments_%g_a_%g_kernel_size_%g_%g.eps', nor, calib_size(1), calib_size(2), ssg.section, ssg.alpha, ssg.patch_size(1), ssg.patch_size(2)))

end


if para.setting.ifplot && para.Recon.ssg.section == 8
    ssg = para.Recon.ssg;
    
    nor = 300;
    N_test = NUFFT.init(Data.kx(:,1:nor,:), Data.ky(:,1:nor,:), 1, [6,6], sx, sx);
    img_test = NUFFT.NUFFT_adj(Data.kSpace(:,1:nor,:,:).*conj(Data.phase_mod(:,1:nor,:,:,:)), N_test);
    img_test_sens = sum(img_test .* conj(Data.sens_map), 4);
    img_test = section_kernel_apply_8_seg_2(img_test, Data.ssg);
    img_test_ssg_sens = squeeze(sum(img_test.*conj(Data.sens_map), 4));
    
    figure
    imshow = abs([fliplr(img_test_sens(:,:,1,:)); fliplr(img_test_ssg_sens(:,:,1,:)); fliplr(abs(img_test_sens(:,:,1,:)- img_test_ssg_sens(:,:,1,:)))]);
    imshow = imshow(:,:);
    imagesc(imshow);axis image
    colormap gray; brighten(0.4); axis off
    title(sprintf('%g ray/frame, calibration size = [%g,%g],\n k-space segments = %g, \\alpha = %g, kernel size = [%g,%g]', nor, calib_size(1), calib_size(2), ssg.section, ssg.alpha, ssg.patch_size(1), ssg.patch_size(2)), 'FontSize', 30)
    set(gcf, 'Position', [0,0,900,900]);

end



%% PSF
if para.setting.ifplot
    ssg = para.Recon.ssg;
    
    nor = 300;
    psf_position_x = sx/2;
    psf_position_y = sx/2;
    Colors = colors;
    
    img_psf = zeros(sx, sx, nof_total);
    img_psf(psf_position_x, psf_position_y, :) = 1;
   
   
    N_test = NUFFT.init(Data.kx(:,1:nor,:), Data.ky(:,1:nor,:), 1, [6,6], sx, sx);
    
    kspace = NUFFT.NUFFT(img_psf.*Data.sens_map, N_test);
    kspace = kspace.*Data.phase_mod(:,1:nor,:,:,:);
    kspace = sum(kspace, 5);
    
    img_test = NUFFT.NUFFT_adj(kspace .* conj(Data.phase_mod(:,1:nor,:,:,:)), N_test);
    
    img_test_sens = sum(img_test .* conj(Data.sens_map),4);
    
    img_test_ssg_sens = section_kernel_apply(img_test, Data.ssg);
    img_test_ssg_sens = sum(img_test_ssg_sens .* conj(Data.sens_map), 4);
    
    f = figure;
    imshow = abs([fliplr(img_test_sens(:,:,1,:)); fliplr(img_test_ssg_sens(:,:,1,:)); fliplr(abs(img_test_sens(:,:,1,:)- img_test_ssg_sens(:,:,1,:)))]);
    imshow = imshow(:,:);
    imagesc(imshow);axis image
    colormap gray; brighten(0.6); axis off
    title(sprintf('PSF: %g ray/frame, calibration size = [%g,%g],\n k-space segments = %g, \\alpha = %g, kernel size = [%g,%g]', nor, calib_size(1), calib_size(2), ssg.section, ssg.alpha, ssg.patch_size(1), ssg.patch_size(2)), 'FontSize', 30)
    set(gcf, 'Position', [0,0,900,900])
    set(gca, 'pos', [0.05 0.05 0.9 0.85])
    text(50, 10, 'w/o SSG', 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center', 'Color', 'white', 'FontSize', 20)
    text(50, 10+sx, 'w/ SSG', 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center', 'Color', 'white', 'FontSize', 20)
    text(50, 10+sx*2, 'difference', 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center', 'Color', 'white', 'FontSize', 20)
    text(150, 10, 'slice 1', 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center', 'Color', colors(1), 'FontSize', 20)
    text(150+sx, 10, 'slice 2', 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center', 'Color', colors(2), 'FontSize', 20)
    text(150+sx*2, 10, 'slice 3', 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center', 'Color', colors(3), 'FontSize', 20)
    patch([1,sx-2,sx-2,1],[0,0,sx*3+1,sx*3+1],Colors(1,:),'FaceColor','none', 'linewidth', 4, 'Edgecolor', Colors(1,:))
    patch([sx+2,sx*2-2,sx*2-2,sx+2],[0,0,sx*3+1,sx*3+1],Colors(1,:),'FaceColor','none', 'linewidth', 4, 'Edgecolor', Colors(2,:))
    patch([sx*2+2,sx*3,sx*3,sx*2+2],[0,0,sx*3+1,sx*3+1],Colors(1,:),'FaceColor','none', 'linewidth', 4, 'Edgecolor', Colors(3,:))
    filename = sprintf('./Figure/PSF_%g_ray_calibration_size_%g_%g_segments_%g_alpha_%g_kernel_size_%g_%g.eps', nor, calib_size(1), calib_size(2), ssg.section, ssg.alpha, ssg.patch_size(1), ssg.patch_size(2));
    hgexport(f, filename)
    
    f = figure;
    hold on
    line_x_sens = squeeze(img_test_sens(psf_position_x,:,:,:,:,:));
    line_x_ssg_sens = squeeze(img_test_ssg_sens(psf_position_x,:,:,:,:,:));
    for i=1:nsms
        plot(-sx/2:sx/2-1,squeeze(abs(line_x_ssg_sens(:,1,i))), 'LineWidth', 2, 'Color', Colors(i,:))
        plot(-sx/2:sx/2-1,squeeze(abs(line_x_sens(:,1,i))), 'LineWidth', 2, 'LineStyle', '--', 'Color', Colors(i,:))
    end
    title(sprintf('PSF (y = 0): %g ray/frame, calibration size = [%g,%g],\n k-space segments = %g, \\alpha = %g, kernel size = [%g,%g]', nor, calib_size(1), calib_size(2), ssg.section, ssg.alpha, ssg.patch_size(1), ssg.patch_size(2)), 'FontSize', 30)
    axis([-sx/2,sx/2-1,0,0.01])
    set(gcf, 'Position', [0,0,900,400]);
    set(gca,'pos',[0.1 0.1 0.85 0.7])
    set(gca, 'FontSize', 26)
    box on
    grid on
    legend('slice 1 w/ SSG', 'slice 1 w/o SSG', 'slice 2 w/ SSG', 'slice 2 w/o SSG', 'slice 3 w/ SSG', 'slice 3 w/o SSG')
    filename = sprintf('./Figure/PSF_x_%g_ray_calibration_size_%g_%g_segments_%g_alpha_%g_kernel_size_%g_%g.eps', nor, calib_size(1), calib_size(2), ssg.section, ssg.alpha, ssg.patch_size(1), ssg.patch_size(2));
    hgexport(f, filename)
    
    f = figure;
    hold on
    line_y_sens = squeeze(img_test_sens(:,psf_position_y,:,:,:,:));
    line_y_ssg_sens = squeeze(img_test_ssg_sens(:,psf_position_y,:,:,:,:));
    for i=1:nsms
        plot(-sx/2:sx/2-1,squeeze(abs(line_y_ssg_sens(:,1,i))), 'LineWidth', 2, 'Color', Colors(i,:))
        plot(-sx/2:sx/2-1,squeeze(abs(line_y_sens(:,1,i))), 'LineWidth', 2, 'LineStyle', '--', 'Color', Colors(i,:))
    end
    title(sprintf('PSF (x = 0): %g ray/frame, calibration size = [%g,%g],\n k-space segments = %g, \\alpha = %g, kernel size = [%g,%g]', nor, calib_size(1), calib_size(2), ssg.section, ssg.alpha, ssg.patch_size(1), ssg.patch_size(2)), 'FontSize', 30)
    axis([-sx/2,sx/2-1,0,0.01])
    set(gcf, 'Position', [0,0,900,400]);
    set(gca,'pos',[0.1 0.1 0.85 0.7])
    set(gca, 'FontSize', 26)
    box on
    grid on
    legend('slice 1 w/ SSG', 'slice 1 w/o SSG', 'slice 2 w/ SSG', 'slice 2 w/o SSG', 'slice 3 w/ SSG', 'slice 3 w/o SSG')
    filename = sprintf('./Figure/PSF_y_%g_ray_calibration_size_%g_%g_segments_%g_alpha_%g_kernel_size_%g_%g.eps', nor, calib_size(1), calib_size(2), ssg.section, ssg.alpha, ssg.patch_size(1), ssg.patch_size(2));
    hgexport(f, filename)
    %% PSF slice aliasing
%     nor = 300;
   
%     N_test = NUFFT.init(Data.kx(:,1:nor,:), Data.ky(:,1:nor,:), 1, [6,6], sx, sx);
    
    kspace = NUFFT.NUFFT(img_psf.*Data.sens_map, N_test);
    kspace = kspace.*Data.phase_mod(:,1:nor,:,:,:);
%     kspace = sum(kspace, 5);
    
    img_test = NUFFT.NUFFT_adj(kspace, N_test);
    
    img_test_sens = sum(img_test .* conj(Data.sens_map),4);
    
    img_test_ssg_sens = section_kernel_apply(img_test, Data.ssg);
    img_test_ssg_sens = sum(img_test_ssg_sens .* conj(Data.sens_map), 4);
    
    f = figure;
    imshow = abs([fliplr(img_test_sens(:,:,1,:)); fliplr(img_test_ssg_sens(:,:,1,:)); fliplr(abs(img_test_sens(:,:,1,:)- img_test_ssg_sens(:,:,1,:)))]);
    imshow = imshow(:,:);
    imagesc(imshow);axis image
    colormap gray; brighten(0.6); axis off
    title(sprintf('PSF %g ray/frame, calibration size = [%g,%g],\n k-space segments = %g, \\alpha = %g, kernel size = [%g,%g]', nor, calib_size(1), calib_size(2), ssg.section, ssg.alpha, ssg.patch_size(1), ssg.patch_size(2)), 'FontSize', 30)
    set(gcf, 'Position', [0,0,900,900])
    set(gca, 'pos', [0.05 0.05 0.9 0.85])
    text(50, 10, 'w/o SSG', 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center', 'Color', 'white', 'FontSize', 20)
    text(50, 10+sx, 'w/ SSG', 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center', 'Color', 'white', 'FontSize', 20)
    text(50, 10+sx*2, 'difference', 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center', 'Color', 'white', 'FontSize', 20)
    text(150, 10, 'slice 1', 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center', 'Color', colors(1), 'FontSize', 20)
    text(100+sx, 10, 'leakage to slice 2', 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center', 'Color', colors(2), 'FontSize', 20)
    text(100+sx*2, 10, 'leakage to slice 3', 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center', 'Color', colors(3), 'FontSize', 20)
    patch([1,sx-2,sx-2,1],[0,0,sx*3+1,sx*3+1],Colors(1,:),'FaceColor','none', 'linewidth', 4, 'Edgecolor', Colors(1,:))
    patch([sx+2,sx*2-2,sx*2-2,sx+2],[0,0,sx*3+1,sx*3+1],Colors(1,:),'FaceColor','none', 'linewidth', 4, 'Edgecolor', Colors(2,:))
    patch([sx*2+2,sx*3,sx*3,sx*2+2],[0,0,sx*3+1,sx*3+1],Colors(1,:),'FaceColor','none', 'linewidth', 4, 'Edgecolor', Colors(3,:))
    filename = sprintf('./Figure/PSF_leakage_%g_ray_calibration_size_%g_%g_segments_%g_alpha_%g_kernel_size_%g_%g.eps', nor, calib_size(1), calib_size(2), ssg.section, ssg.alpha, ssg.patch_size(1), ssg.patch_size(2));
    hgexport(f, filename)
    
    f = figure;
    hold on
    line_x_sens = squeeze(img_test_sens(sx/2,:,:,:,:,:));
    line_x_ssg_sens = squeeze(img_test_ssg_sens(sx/2,:,:,:,:,:));
    for i=1:nsms
        plot(-sx/2:sx/2-1,squeeze(abs(line_x_ssg_sens(:,1,i))), 'LineWidth', 2, 'Color', Colors(i,:))
        plot(-sx/2:sx/2-1,squeeze(abs(line_x_sens(:,1,i))), 'LineWidth', 2, 'LineStyle', '--', 'Color', Colors(i,:))
    end
    title(sprintf('PSF (y = 0): %g ray/frame, calibration size = [%g,%g],\n k-space segments = %g, \\alpha = %g, kernel size = [%g,%g]', nor, calib_size(1), calib_size(2), ssg.section, ssg.alpha, ssg.patch_size(1), ssg.patch_size(2)), 'FontSize', 30)
    axis([-sx/2,sx/2-1,0,0.01])
    set(gcf, 'Position', [0,0,900,400]);
    set(gca,'pos',[0.1 0.1 0.85 0.7])
    set(gca, 'FontSize', 26)
    box on
    grid on
    legend('target slice w/ SSG', 'target slice w/o SSG', 'aliased slice 1 w/ SSG', 'aliased slice 1 w/o SSG', 'aliased slice 2 w/ SSG', 'aliased slice 2 w/o SSG')
    filename = sprintf('./Figure/PSF_leakage_x_%g_ray_calibration_size_%g_%g_segments_%g_alpha_%g_kernel_size_%g_%g.eps', nor, calib_size(1), calib_size(2), ssg.section, ssg.alpha, ssg.patch_size(1), ssg.patch_size(2));
    hgexport(f, filename)
    
    
    f = figure;
    hold on
    line_y_sens = squeeze(img_test_sens(:,sx/2,:,:,:,:));
    line_y_ssg_sens = squeeze(img_test_ssg_sens(:,sx/2,:,:,:,:));
    for i=1:nsms
        plot(-sx/2:sx/2-1,squeeze(abs(line_y_ssg_sens(:,1,i))), 'LineWidth', 2, 'Color', Colors(i,:))
        plot(-sx/2:sx/2-1,squeeze(abs(line_y_sens(:,1,i))), 'LineWidth', 2, 'LineStyle', '--', 'Color', Colors(i,:))
    end
    title(sprintf('PSF (x = 0): %g ray/frame, calibration size = [%g,%g],\n k-space segments = %g, \\alpha = %g, kernel size = [%g,%g]', nor, calib_size(1), calib_size(2), ssg.section, ssg.alpha, ssg.patch_size(1), ssg.patch_size(2)), 'FontSize', 30)
    axis([-sx/2,sx/2-1,0,0.01])
    set(gcf, 'Position', [0,0,900,400]);
    set(gca,'pos',[0.1 0.1 0.85 0.7])
    set(gca, 'FontSize', 26)
    box on
    grid on
    legend('target slice w/ SSG', 'target slice w/o SSG', 'aliased slice 1 w/ SSG', 'aliased slice 1 w/o SSG', 'aliased slice 2 w/ SSG', 'aliased slice 2 w/o SSG')
    filename = sprintf('./Figure/PSF_leakage_y_%g_ray_calibration_size_%g_%g_segments_%g_alpha_%g_kernel_size_%g_%g.eps', nor, calib_size(1), calib_size(2), ssg.section, ssg.alpha, ssg.patch_size(1), ssg.patch_size(2));
    hgexport(f, filename)
end

end


% [kernel, kernel_im] = SSG_train_section_main(kSpace_calib_source, para);

% sens_conj_ssg = squeeze(sum(permute(conj(Data.sens_map),[1,2,3,6,4,5]).*kernel_im,5));
% sens_conj_ssg = permute(sens_conj_ssg,[1,2,5,3,4]);