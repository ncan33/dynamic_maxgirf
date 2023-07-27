function [sweep, tTV_grid, sTV_grid] = dual_te_STCR_parameter_sweep(narm_frame, tTV_step_factor, sTV_step_factor, ...
    n_tTV_steps, n_sTV_steps, max_sTV, perform_initial_tTV_sweep, tTV_anchor, tTV_low, tTV_high, niter, ...
    ifsave, ifvideo, ifGPU, path)
    % Edit that needs to be made: number of sweeps in temporal direction
    % needs to be a customizable parameter
    % 
    % Parameter sweep for spatiotemporally constrained reconstruction on
    % dual TE variable density spiral raw RTHawk data. It calls STCR on my
    % the data of this repository
    
    arguments
        narm_frame
        tTV_step_factor = 5 % each step will be x5 farther from the anchor
        sTV_step_factor = 5
        n_tTV_steps = 4
        n_sTV_steps = 5
        max_sTV = 0.05
        perform_initial_tTV_sweep = 0
        tTV_anchor = 0.02
        tTV_low = 1e-7
        tTV_high = 1e-1
        niter = 75
        ifsave = 1
        ifvideo = 1
        ifGPU = 1
        path = '/server/sdata/ncan/mri_data/disc/lung/vol0457_20221021/raw_hawk/usc_disc_yt_2022_10_21_133643_dual-te_dynamic.mat'
    end
    
    if n_tTV_steps > 6
        error(['Number of steps in the temporal direction too high.', ...
            ' Ensure the number of steps is less than 6'])
    end
    
    if n_tTV_steps < 4
        error(['Number of steps in the temporal direction too low.', ...
            ' Ensure the number of steps is more than 3'])
    end
    
    if ~(n_tTV_steps >= 4 && n_tTV_steps <= 5)
        input('Only 4 or 5 steps in the temporal direction are recommended. To proceed, press enter.')
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
    
    %% make the initial tTV sweep vector
    n_tTV_iter = log10(tTV_high/tTV_low) + 1;
    initial_tTV_sweep = zeros(1, n_tTV_iter);
    
    for i = 1:n_tTV_iter
        initial_tTV_sweep(i) = tTV_low*10^(i-1);
    end
    
    %% sweep through initial tTVs
    initial_sTV_sweep = 1e-5;
    
    if perform_initial_tTV_sweep == 1
        for i = 1:n_tTV_iter
            [im_echo_1, im_echo_2, NUFFT_im_echo_1, NUFFT_im_echo_2, kspace_info, para] = dual_te_STCR_wrapper(narm_frame, initial_tTV_sweep(i), initial_sTV_sweep, niter, 0, ifGPU, 0);
            if ifsave
                save_name = sprintf(['./recon_data/parameter_sweep/', num2str(narm_frame), 'arm_', num2str(para.weight_tTV), '_tTV_', num2str(para.weight_sTV),'_sTV_','%s_recon.mat'], dir(path).name(1:end-8));
                save(save_name, 'im_echo_1', 'im_echo_2', 'NUFFT_im_echo_1', 'NUFFT_im_echo_2', 'kspace_info', 'para', '-v7.3');
            end
            if i == 1
                initial_sweep = im_echo_1;
            else
                initial_sweep = [initial_sweep, im_echo_1];
            end
        end
    
        close all
        disp('Initial sweep done! Yay!')
        input('Hit enter, and the initial parameter sweep will play.');

        play_mri_video(29, 20/narm_frame, initial_sweep)
        tTV_anchor = input('Enter the best tTV: ');
        clear initial_sweep
    elseif perform_initial_tTV_sweep == 0
        disp(['Using ', num2str(tTV_anchor), ' as the tTV_anchor...'])
    else
        error('perform_initial_tTV_sweep can either be set to true or false')
    end
    
    %% sweep though anchor tTV
    
    [tTV_sweep, sTV_sweep] = generate_anchored_sweep_vectors(n_tTV_steps, ...
    n_sTV_steps, tTV_step_factor, sTV_step_factor, tTV_anchor, max_sTV);
    
    for i = 1:length(tTV_sweep)
        for j = 1:length(sTV_sweep)
            save_name = sprintf(['./recon_data/parameter_sweep/', num2str(narm_frame), 'arm_', num2str(tTV_sweep(i)), '_tTV_', num2str(sTV_sweep(j)),'_sTV_','%s_recon.mat'], dir(path).name(1:end-8));
            if ~isfile(save_name) && ~(tTV_sweep(i)==0 && sTV_sweep(j)==0) % don't STCR again if it already exists and also don't STCR for [tTV,sTV] = [0,0]
                [im_echo_1, im_echo_2, NUFFT_im_echo_1, NUFFT_im_echo_2, kspace_info, para] = dual_te_STCR_wrapper(narm_frame, tTV_sweep(i), sTV_sweep(j), niter, 0, ifGPU, 0);
                save(save_name, 'im_echo_1', 'NUFFT_im_echo_1', 'kspace_info', 'para', '-v7.3');
                disp('Save complete')
            else
                disp('File already exists')
            end
        end
    end
    
    for i = 1:length(tTV_sweep)
        disp('Achored sweep done! This is great!')
    end
    
    %% imtile the data -- echo 1 only for simplicity
    clearvars -except narm_frame tTV_step_factor sTV_step_factor tTV_anchor n_tTV_steps n_sTV_steps ifsave ifvideo path
    
    [sweep, tTV_grid, sTV_grid] = parameter_sweep_read_only(narm_frame, tTV_step_factor, ...
        sTV_step_factor, tTV_anchor, n_tTV_steps, n_sTV_steps, ifsave, ifvideo, path);
    
end