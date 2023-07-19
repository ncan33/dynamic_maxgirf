function sweep = dual_te_STCR_parameter_sweep(narm_frame, tTV_low, tTV_high, ifGPU, path)
    % Parameter sweep for spatiotemporally constrained reconstruction on
    % dual TE variable density spiral raw RTHawk data. It calls STCR on my
    % the data of this repository
    
    arguments
        narm_frame
        tTV_low = 0.0001
        tTV_high = 0.1
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
    
    %% make a vector
    n_tTV_iter = log10(tTV_high/tTV_low) + 1;
    initial_tTV_sweep = zeros(1, n_tTV_iter);
    
    for i = 1:n_tTV_iter
        initial_tTV_sweep(i) = tTV_low*10^(n_tTV_iter-1);
    end
    disp(initial_tTV_sweep)
    
    for i = 1:n_tTV_iter
        %disp(i)
        %[im_echo_1, im_echo_2, NUFFT_im_echo_1, NUFFT_im_echo_2, kspace_info, para] = dual_te_STCR_wrapper(narm_frame,tTV,sTV,ifGPU,path)
    end
end
    