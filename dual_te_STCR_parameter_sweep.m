function sweep = dual_te_STCR_parameter_sweep(narm_frame, tTV_step_factor, sTV_step_factor, ...
    n_tTV_steps, n_sTV_steps, tTV_low, tTV_high, niter, ifsave, ifGPU, path)
    % Edit that needs to be made: number of sweeps in temporal direction
    % needs to be a customizable parameter
    % 
    % Parameter sweep for spatiotemporally constrained reconstruction on
    % dual TE variable density spiral raw RTHawk data. It calls STCR on my
    % the data of this repository
    
    arguments
        narm_frame
        tTV_step_factor = 10 % each step will be x10 farther from the anchor
        sTV_step_factor = 10
        n_tTV_steps = 4
        tTV_low = 1e-7
        tTV_high = 1e-1
        niter = 75
        %niter = 150 % this is plenty
        ifsave = 1
        ifGPU = 1
        path = '/server/sdata/ncan/mri_data/disc/lung/vol0457_20221021/raw_hawk/usc_disc_yt_2022_10_21_133643_dual-te_dynamic.mat'
    end
    
    if n_tTV_steps > 6
        error(['Number of steps in the temporal direction too high.' ...
            ' Ensure the number of steps is below'])
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
    %input('Hit enter, and the initial parameter sweep will play.');
    %UNCOMMENT THE ABOVE
    %UNVOMMENT THE ABOVE
    play_mri_video(29, 20/narm_frame, initial_sweep)
    %tTV_anchor = input('Enter the best tTV: ');
    %UNCOMMENT THE ABOVE
    %UNVOMMENT THE ABOVE
    %comment the below
    tTV_anchor = 1e-1;
    clear initial_sweep
    
    %% sweep though anchor tTV
    
    tTV_sweep = zeros(1,n_);
    tTV_sweep(3) = tTV_anchor; tTV_sweep(1) = tTV_anchor/(tTV_step_factor^2);
    tTV_sweep(2) = tTV_anchor/tTV_step_factor; tTV_sweep(4) = tTV_anchor*tTV_step_factor;
    
    sTV_sweep = zeros(1,n_sTV_steps);
    for i = 1:n_sTV_steps
        sTV_sweep(i) = 1*sTV_step_factor^(-n_sTV_steps)*sTV_step_factor^(i-1);
    end
    
    for i = 1:length(tTV_sweep)
        for j = 1:length(sTV_sweep)
            [im_echo_1, im_echo_2, NUFFT_im_echo_1, NUFFT_im_echo_2, kspace_info, para] = dual_te_STCR_wrapper(narm_frame, tTV_sweep(i), initial_sTV_sweep, niter, 0, ifGPU, 0);
            if ifsave
                save_name = sprintf(['./recon_data/parameter_sweep/', num2str(narm_frame), 'arm_', num2str(tTV_sweep(i)), '_tTV_', num2str(sTV_sweep(j)),'_sTV_','%s_recon.mat'], dir(path).name(1:end-8));
                save(save_name, 'im_echo_1', 'im_echo_2', 'NUFFT_im_echo_1', 'NUFFT_im_echo_2', 'kspace_info', 'para', '-v7.3');
            end
        end
    end
    
    for i = 1:length(tTV_sweep)
        disp('Achored sweep done! This is great!')
    end
    
    %% imtile the data -- echo 1 only for simplicity
    clearvars -except tTV_sweep sTV_sweep narm_frame ifsave path
    
    for i = 1:length(tTV_sweep)
        for j = 1:length(sTV_sweep)
            load_name = sprintf(['./recon_data/parameter_sweep/', num2str(narm_frame), 'arm_', num2str(tTV_sweep(i)), '_tTV_', num2str(sTV_sweep(j)),'_sTV_','%s_recon.mat'], dir(path).name(1:end-8));
            load(load_name, 'im_echo_1');
            
            if i == 1 && j == 1
                sweep_row_1 = im_echo_1;
            elseif i == 1
                sweep_row_1 = [sweep_row_1, im_echo_1];
            elseif i == 2 && j == 1
                disp('Successfully sweeped through first row!')
                sweep_row_2 = im_echo_1;
            elseif i == 2
                sweep_row_2 = [sweep_row_2, im_echo_1];
            elseif i == 3 && j == 1
                disp('Successfully sweeped through second row!')
                sweep_row_3 = im_echo_1;
            elseif i == 3
                sweep_row_3 = [sweep_row_3, im_echo_1];
            elseif i == 4 && j == 1
                disp('Successfully sweeped through third row!')
                sweep_row_4 = im_echo_1;
            else
                sweep_row_4 = [sweep_row_4, im_echo_1];
            end
        end
    end
    disp('Successfully sweeped through fourth row!')
    
    sweep = [sweep_row_1; sweep_row_2; sweep_row_3; sweep_row_4];
    if ifsave
        save(['sweep_',num2str(narm_frame),'_arm'],'sweep')
    end
end
    