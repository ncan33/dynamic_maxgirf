function [PSF, PSF_cropped] = RTHawk_spiral_PSF_visualizer(narm_frame, kspace_info, fov_factor, ifsave)
    % Generate PSF of the spiral trajectory for RTHawk data for a user
    % defined number of arms per frame. Uses NUFFT.
    % 
    % Note!
    % This function is for dual TE data, and it performs splitting on the
    % kspace trajectory. If single TE data exists, splitting is not
    % necessary. A future update to this function will add dual TE as a
    % 'case' rather than a built-in functionality.
    arguments
        narm_frame
        kspace_info = []
        fov_factor = 1
        ifsave = 0
    end
    
    if isempty(kspace_info)
        load('supplementary/kspace_info.mat', 'kspace_info');
    else
        disp('Using custom kspace_info input by user...')
    end
    %% add paths
    run('dynamic_maxgirf_setup.m')
    
    %% assign kspace attributes to variables
    nview = kspace_info.user_interleaves; % total number of interleaves in original kspace trajectory
    
    view_order = kspace_info.viewOrder; % view order

    echo_idx = view_order > nview;
    echo_idx = echo_idx + 1; % if the element in echo_idx == 1, then echo_1,
                             % and when echo_idx == 2, then echo_2
                             
    nsample = kspace_info.kspace.samples; % number of samples
    ncoil = kspace_info.extent(2); % number of coils
    
    res = [kspace_info.user_ResolutionX, kspace_info.user_ResolutionY];
    fov = [kspace_info.user_FieldOfViewX, kspace_info.user_FieldOfViewY];
    
    matrix_size_keep = round(fov ./ res / 2) * 2;
    matrix_size = round(matrix_size_keep * 2 / 2) * fov_factor;
    para.Recon.matrix_size = matrix_size;
    para.Recon.FOV = fov(1)/100; % units of decimeter for some reason
    
    %% organize kspace trajectory
    narm_total = floor(length(kspace_info.viewOrder) / 2); % number of arms per echo
    nframes = floor(narm_total / narm_frame);
    narm_total = nframes * narm_frame;
    
    view_order = view_order(echo_idx == 1); % view_order for echo_1
    view_order(narm_total + 1 : end) = []; % discard excess views
    view_order  = reshape(view_order, [narm_frame, nframes]);
    
    kx = kspace_info.kx_GIRF; % kx
    ky = kspace_info.ky_GIRF; % ky
    
    kx_echo_1 = zeros(nsample, narm_frame, nframes); %kx_echo_1 = zeros(size(kspace_echo_1,1), size(kspace_echo_1,2), size(kspace_echo_1,3));
    ky_echo_1 = kx_echo_1;
    
    for i = 1:narm_frame
        for j = 1:nframes
            kx_echo_1(:,i,j) = kx(:, view_order(i,j));
            ky_echo_1(:,i,j) = ky(:, view_order(i,j));
        end
    end
    
    kx = kx_echo_1;
    ky = ky_echo_1;
    
    %% perform NUFFT
    Data.N = NUFFT.init(kx*para.Recon.matrix_size(1), ky*para.Recon.matrix_size(2), 1, [4, 4], para.Recon.matrix_size(1), para.Recon.matrix_size(1));
    Data.N.W = kspace_info.DCF(:, 1);
    
    Data.kSpace = ones([size(kx), ncoil]);
    Data.first_est = NUFFT.NUFFT_adj(Data.kSpace, Data.N);
    
    NUFFT_im = sum(bsxfun(@times, conj(get_sens_map(Data.first_est, '2D')), Data.first_est), 4);
    NUFFT_im = fliplr(rot90(NUFFT_im, -1)); % adjust orientation
    
    PSF = NUFFT_im;
    clear NUFFT_im
    
    disp('NUFFT complete.')
    
    %% normalize the PSF
    last_frame = PSF(:, :, end);
    sz = size(last_frame, 1);
    center_pixel_mag = abs(last_frame(sz/2+1, sz/2+1));
    PSF = PSF / center_pixel_mag;
    
    %% plot PSF image and save it
    
    % declare axes
    axis_data_high = (para.Recon.FOV * 10 / 2) * fov_factor; % upper bound for new XTick vector
    axis_data_low = - axis_data_high; % lower bound for new XTick vector
    axis_data = round(linspace(axis_data_low, axis_data_high, 200)); % new XTickLabel vector
    % plot
    imagesc('XData', axis_data, 'YData', axis_data, 'CData', log10(abs(PSF(:,:,end))))
    caxis([-5.5 0])
    xlabel('X (cm)'); ylabel('Y (cm)'); axis image; colorbar; colormap gray
    title(['Log-scale PSF for ', num2str(narm_frame), ' arms per frame'])
    if ifsave
        saveas(gcf, ['./figures/', num2str(narm_frame), 'arm_logscale_PSF.png'])
    end
    
    % zoom into mainlobe
    sz = size(PSF, 1);
    zoom_amount = 10; % number of pixels
    PSF_cropped = PSF((sz / 2)-zoom_amount+1:(sz / 2)+zoom_amount+1, ...
        (sz / 2)-zoom_amount+1:(sz / 2)+zoom_amount+1, :);
    % declare axes
    axis_data_high = (para.Recon.FOV * 100 / 2) * fov_factor * size(PSF_cropped, 1) / sz; % upper bound for new XTick vector
    axis_data_low = - axis_data_high; % lower bound for new XTick vector
    axis_data = round(linspace(axis_data_low, axis_data_high, 50)); % new XTickLabel vector
    % plot
    imagesc(axis_data, axis_data, log10(abs(PSF_cropped(:,:,end))));
    caxis([-2 0])
    xlabel('X (mm)'); ylabel('Y (mm)'); axis image; colorbar; colormap gray
    title(['Mainlobe of PSF for ', num2str(narm_frame), ' arms per frame (log-scale)'])
    if ifsave
        saveas(gcf, ['./figures/', num2str(narm_frame), 'arm_mainlobe_of_PSF.png'])
    end
    
end