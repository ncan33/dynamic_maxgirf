function sweep = parameter_sweep_read_only(narm_frame, tTV_step_factor, sTV_step_factor, tTV_anchor, tTV_high, ...
    niter, ifsave, ifGPU, path)
    
    % FUNCTION MUST BE RENAMED AND RELOCATED
    % 
    % Read-only parameter sweep for spatiotemporally constrained
    % reconstruction on dual TE variable density spiral raw RTHawk data.
    % It calls STCR on my the data of this repository.
    
    arguments
        narm_frame
        tTV_step_factor = 10 % each step will be x10 farther from the anchor
        sTV_step_factor = 10
        tTV_anchor = 1e-3
        ifsave = 1
        path = '/server/sdata/ncan/mri_data/disc/lung/vol0457_20221021/raw_hawk/usc_disc_yt_2022_10_21_133643_dual-te_dynamic.mat'
    end
    
    %% add paths
    addpath ./util/mfile/functions/
    addpath ./util/mfile/registrtation/
    addpath ./util/mfile/quantification/
    addpath ./util/mfile/vdspiral/
    addpath ./util/
    
    if ~isfolder('./recon_data/parameter_sweep')
        error("You don't have any data to read")
    end
    
    %% create sweep vectors
    tTV_sweep = zeros(1,4);
    tTV_sweep(3) = tTV_anchor; tTV_sweep(1) = tTV_anchor/(tTV_step_factor^2);
    tTV_sweep(2) = tTV_anchor/tTV_step_factor; tTV_sweep(4) = tTV_anchor*tTV_step_factor;

    sTV_sweep = zeros(1,6);
    for i = 1:6
        sTV_sweep(i) = 1e-6*sTV_step_factor^(i-1);
    end
    
    %% sweep through data
    for i = 1:length(tTV_sweep)
        for j = 1:length(sTV_sweep)
            load_name = sprintf(['./recon_data/parameter_sweep/', num2str(narm_frame), 'arm_', num2str(tTV_sweep(i)), '_tTV_', num2str(sTV_sweep(j)),'_sTV_','%s_recon.mat'], dir(path).name(1:end-8));
            load(load_name, 'im_echo_1');
            
            if i == 1 && j == 1
                sweep_row_1 = im_echo_1;
            elseif i == 1
                sweep_row_1 = [sweep_row_1, im_echo_1];
            elseif i == 2 && j == 1
                sweep_row_2 = im_echo_1;
            elseif i == 2
                sweep_row_2 = [sweep_row_2, im_echo_1];
            elseif i == 3 && j == 1
                sweep_row_3 = im_echo_1;
            elseif i == 3
                sweep_row_3 = [sweep_row_3, im_echo_1];
            elseif i == 4 && j == 1
                sweep_row_4 = im_echo_1;
            else
                sweep_row_4 = [sweep_row_4, im_echo_1];
            end
        end
    end
    
    disp(size(sweep_row_1))
    disp(size(sweep_row_2))
    disp(size(sweep_row_3))
    disp(size(sweep_row_4))
    
    sweep = [sweep_row_1; sweep_row_2; sweep_row_3, sweep_row_4];
    if ifsave
        save(['sweep_',num2str(narm_frame),'_arm'],'sweep')
    end
end