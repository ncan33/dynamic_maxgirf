function [im_echo, NUFFT_im, para] = STCR(kspace_info, kspace, kx, ky, para)
    %% matrix size
    res = [kspace_info.user_ResolutionX, kspace_info.user_ResolutionY];
    fov = [kspace_info.user_FieldOfViewX, kspace_info.user_FieldOfViewY];
    matrix_size_keep    = round(fov ./ res / 2) * 2;
    matrix_size         = round(matrix_size_keep * 2 / 2) * 1.5;
    para.Recon.matrix_size = matrix_size;

    %% perform STCR for echo 1
    if isnumeric(para.Recon.time_frames)
        time_frames = para.Recon.time_frames;

        kx = kx(:, :, time_frames);
        ky = ky(:, :, time_frames);
        kspace = kspace(:, :, time_frames, :);
    end

    para.Recon.FOV = fov(1)/100; %units of decimeter for some reason
    
    Data.N = NUFFT.init(kx*para.Recon.matrix_size(1), ky*para.Recon.matrix_size(2), 1, [4, 4], para.Recon.matrix_size(1), para.Recon.matrix_size(1));
    Data.N.W = kspace_info.DCF(:, 1);

    Data.kSpace = kspace;
    Data.first_est = NUFFT.NUFFT_adj(Data.kSpace, Data.N);
    
    NUFFT_im = sum(bsxfun(@times, conj(get_sens_map(Data.first_est, '2D')), Data.first_est), 4);
    
    disp('NUFFT complete.')
    
    scale = max(abs(Data.first_est(:)));

    Data.sens_map = get_sens_map(Data.first_est, '2D');
    Data.first_est = sum(Data.first_est .* conj(Data.sens_map), 4);

    para.Recon.weight_tTV = scale * para.weight_tTV; % temporal regularization weight
    para.Recon.weight_sTV = scale * para.weight_sTV; % spatial regularization weight

    disp(['Beginning STCR for tTV = ', num2str(para.weight_tTV), ' and sTV = ', num2str(para.weight_sTV)])
    toc
    [Image_recon, para] = STCR_conjugate_gradient(Data, para);
    Image_recon = fliplr(rot90(Image_recon, -1));
    im_echo = crop_half_FOV(Image_recon, matrix_size_keep);
end