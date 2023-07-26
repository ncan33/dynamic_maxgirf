function [sweep, tTV_grid, sTV_grid] = parameter_sweep_read_only(narm_frame, tTV_step_factor, sTV_step_factor, ...
    tTV_anchor, n_tTV_steps, n_sTV_steps, ifsave, path)
    
    % FUNCTION MUST BE RENAMED AND RELOCATED
    % 
    % Read-only parameter sweep for spatiotemporally constrained
    % reconstruction on dual TE variable density spiral raw RTHawk data.
    % It calls STCR on my the data of this repository.
    
    arguments
        narm_frame
        tTV_step_factor = 2.5 % each step will be x10 farther from the anchor
        sTV_step_factor = 10
        tTV_anchor = 1e-1
        n_tTV_steps = 5
        n_sTV_steps = 4
        ifsave = 1
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
        error("You don't have any data to read.")
    end
    
    %% create sweep vectors
    tTV_sweep = zeros(1,n_tTV_steps);
    
    for i = 1:length(tTV_sweep)
        k = -length(tTV_sweep) + 1 + i;
        tTV_sweep(i) = tTV_anchor * (tTV_step_factor^(k));
    end
    
    sTV_sweep = zeros(1,n_sTV_steps);
    for i = 1:n_sTV_steps
        sTV_sweep(i) = 1*sTV_step_factor^(-n_sTV_steps)*sTV_step_factor^(i-1);
    end
    
    tTV_sweep = [0, tTV_sweep]; % add zero column
    sTV_sweep = [0, sTV_sweep]; % add zero column
    
    disp('tTV_sweep is:')
    disp(tTV_sweep)
    disp('sTV_sweep is:')
    disp(sTV_sweep)
    
    %% sweep through data
    tTV_grid = zeros(n_tTV_steps, n_sTV_steps);
    sTV_grid = zeros(n_tTV_steps, n_sTV_steps);
    
    for i = 1:n_tTV_steps % note: n_tTV_steps + 1 == length(tTV_sweep)
        for j = 1:n_sTV_steps  % note: n_sTV_steps + 1 == length(sTV_sweep)
            load_name = sprintf(['./recon_data/parameter_sweep/', num2str(narm_frame), 'arm_', num2str(tTV_sweep(i+1)), '_tTV_', num2str(sTV_sweep(j+1)),'_sTV_','%s_recon.mat'], dir(path).name(1:end-8));
            load(load_name, 'im_echo_1');
           
            tTV_grid(i,j) = tTV_sweep(i+1);
            sTV_grid(i,j) = sTV_sweep(j+1);
            
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
            elseif i == 4
                sweep_row_4 = [sweep_row_4, im_echo_1];
            elseif i == 5 && j == 1
                sweep_row_5 = im_echo_1;
                disp('Successfully sweeped through fourth row!')
            elseif i == 5
                sweep_row_5 = [sweep_row_5, im_echo_1];
            elseif i == 6 && j == 1
                sweep_row_6 = im_echo_1;
                disp('Successfully sweeped through fifth row!')
            elseif i == 6
                sweep_row_6 = [sweep_row_6, im_echo_1];
            end
        end
    end
    disp('Successfully sweeped through last row!')
    
    if n_tTV_steps == 4
        sweep = [sweep_row_1; sweep_row_2; sweep_row_3; sweep_row_4];
    elseif n_tTV_steps == 5
        sweep = [sweep_row_1; sweep_row_2; sweep_row_3; sweep_row_4; sweep_row_5];
    elseif n_tTV_steps == 6
        sweep = [sweep_row_1; sweep_row_2; sweep_row_3; sweep_row_4; sweep_row_5; sweep_row_6];
    end
    
    % completing sweep with the zero column (sTV == 0) and zero row
    % (tTV == 0) and NUFFT    
    for i = 1:n_tTV_steps
        load_name = sprintf(['./recon_data/parameter_sweep/', num2str(narm_frame), 'arm_', num2str(tTV_sweep(i+1)), '_tTV_', num2str(0),'_sTV_','%s_recon.mat'], dir(path).name(1:end-8));
        load(load_name, 'im_echo_1');

        if i == 1
            zero_sTV_column = im_echo_1;
        else
            zero_sTV_column = [zero_sTV_column; im_echo_1];
        end
    end
    
    sweep = [zero_sTV_column, sweep];
    tTV_grid = [transpose(tTV_sweep(end-n_tTV_steps+1:end)), tTV_grid];
    sTV_grid = [zeros(n_tTV_steps,1), sTV_grid];
    
    disp(tTV_grid)
    disp(sTV_grid)
    
    disp('Successfully sweeped through zero sTV column!')
    
    for i = 0:n_sTV_steps
        if i == 0
            load_name = sprintf(['./recon_data/parameter_sweep/', num2str(narm_frame), 'arm_', num2str(0), '_tTV_', num2str(sTV_sweep(i+2)),'_sTV_','%s_recon.mat'], dir(path).name(1:end-8));
            load(load_name, 'im_echo_1');
            load(load_name, 'NUFFT_im_echo_1');
            NUFFT_im_echo_1 = NUFFT_im_echo_1 / mean(max(abs(NUFFT_im_echo_1), [], [1,2]));
            NUFFT_im_echo_1 = NUFFT_im_echo_1 * mean(max(abs(im_echo_1), [], [1,2]));
        else
            load_name = sprintf(['./recon_data/parameter_sweep/', num2str(narm_frame), 'arm_', num2str(0), '_tTV_', num2str(sTV_sweep(i+1)),'_sTV_','%s_recon.mat'], dir(path).name(1:end-8));
            load(load_name, 'im_echo_1');
        end

        if i == 0
            zero_tTV_row = NUFFT_im_echo_1;
        else
            zero_tTV_row = [zero_tTV_row, im_echo_1];
        end
    end
    
    sweep = [zero_tTV_row; sweep];
    tTV_grid = [zeros(1,n_sTV_steps+1); tTV_grid];
    sTV_grid = [sTV_sweep; sTV_grid];
    
    disp('Successfully sweeped through zero tTV column!')
    
    if ifsave
        save(['sweep_',num2str(narm_frame),'_arm'],'sweep','tTV_grid','sTV_grid')
        disp('Successfully saved the sweep variable!')
    end
end