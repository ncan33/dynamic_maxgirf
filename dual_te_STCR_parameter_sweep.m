function sweep = dual_te_STCR_parameter_sweep(narm_frame, step_factor, tTV_low, tTV_high, niter, ifsave, ifGPU, path)
    % Parameter sweep for spatiotemporally constrained reconstruction on
    % dual TE variable density spiral raw RTHawk data. It calls STCR on my
    % the data of this repository
    
    arguments
        narm_frame
        step_factor = 10 % each step will be x10 farther from the anchor
        tTV_low = 1e-6
        tTV_high = 1e-1
        niter = 75
        %niter = 150 % this is plenty
        ifsave = 1
        ifGPU = 1
        path = '/server/sdata/ncan/mri_data/disc/lung/vol0457_20221021/raw_hawk/usc_disc_yt_2022_10_21_133643_dual-te_dynamic.mat'
    end
    
    %% add paths
    addpath ./util/mfile/functions/
    addpath ./util/mfile/registrtation/
    addpath ./util/mfile/quantification/
    addpath ./util/mfile/vdspiral/
    addpath ./util/
    
    if ~isfolder('./recon_data/parameter_sweep')
        mkdir ./recon_data/parameter_sweep
    end
    
    %% make a tTV sweep vector
    n_tTV_iter = log10(tTV_high/tTV_low) + 1;
    initial_tTV_sweep = zeros(1, n_tTV_iter);
    
    for i = 1:n_tTV_iter
        initial_tTV_sweep(i) = tTV_low*10^(i-1);
    end
    
    %% sweep through initial tTVs
    initial_sTV_sweep = 1e-5;
    
    for i = 1:n_tTV_iter
        [im_echo_1, im_echo_2, NUFFT_im_echo_1, NUFFT_im_echo_2, kspace_info, para] = dual_te_STCR_wrapper(narm_frame, initial_tTV_sweep(i), initial_sTV_sweep, niter, 0, ifGPU, 0);
        if ifsave
            save_name = sprintf(['./recon_data/parameter_sweep', num2str(narm_frame), 'arm_', num2str(para.weight_tTV), '_tTV_', num2str(para.weight_sTV),'_sTV_','%s_recon.mat'], all_dat(file_index).name(1:end-8));
            save(save_name, 'im_echo_1', 'im_echo_2', 'NUFFT_im_echo_1', 'NUFFT_im_echo_2', 'kspace_info', 'para', '-v7.3');
        end
    end
    
    tTV_anchor = input('Enter the best tTV: ');
    
    %% sweep again
    
    tTV_sweep = zeros(1,3);
    tTV_sweep(2) = tTV_anchor; tTV_sweep(1) = tTV_anchor/step_factor; tTV_sweep(3) = tTV_anchor*step_factor;
    
    for i = 1:length(tTV_sweep)
        for j = 1:length(sTV_sweep)
            [im_echo_1, im_echo_2, NUFFT_im_echo_1, NUFFT_im_echo_2, kspace_info, para] = dual_te_STCR_wrapper(narm_frame, tTV_sweep(i), initial_sTV_sweep, niter, 0, ifGPU, 0);
            if ifsave
                save_name = sprintf(['./recon_data/parameter_sweep', num2str(narm_frame), 'arm_', num2str(para.weight_tTV), '_tTV_', num2str(para.weight_sTV),'_sTV_','%s_recon.mat'], all_dat(file_index).name(1:end-8));
                save(save_name, 'im_echo_1', 'im_echo_2', 'NUFFT_im_echo_1', 'NUFFT_im_echo_2', 'kspace_info', 'para', '-v7.3');
            end
        end
    end
end
    