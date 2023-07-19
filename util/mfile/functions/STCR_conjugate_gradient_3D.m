function [Image,para] = STCR_conjugate_gradient_3D(Data,para)
% [Image,para] = iteration_recon(first_est,kSpace,sens_map,phase_mod,para,varargin)

disp('Performing iterative STCR reconstruction...');
disp('Showing progress...')
if isempty(Data.first_est)
    Image = [];
    return
end

ifplot         = para.setting.ifplot;
ifGPU          = para.setting.ifGPU;
weight_tTV     = para.Recon.weight_tTV;
weight_sTV     = para.Recon.weight_sTV;
weight_slTV    = para.Recon.weight_sliceTV;
beta_sqrd      = para.Recon.epsilon;

para.Recon.step_size = para.Recon.step_size(1);

new_img_x = single(Data.first_est);

if isfield(Data,'first_guess')
    new_img_x = Data.first_guess;
end

if isfield(Data,'phase_mod')
    Data.phase_mod_conj = conj(single(Data.phase_mod));
end
if isfield(Data,'sens_map')
    Data.sens_map_conj = conj(Data.sens_map);
end

if ifGPU
    new_img_x = gpuArray(new_img_x);
    %Data.kSpace = gpuArray(Data.kSpace);
    
    %Data.first_est = gpuArray(Data.first_est);
    %Data.sens_map = gpuArray(Data.sens_map);
    if isfield(Data,'N')
        for i=1:size(Data.N, 1)
            for j = 1:size(Data.N, 2)
                Data.N(i, j).S = gpuArray(Data.N(i, j).S);
                Data.N(i, j).Apodizer = gpuArray(Data.N(i, j).Apodizer);
                Data.N(i, j).W = gpuArray(Data.N(i, j).W);
            end
        end
    end
    %img_k = gpuArray(img_k);
    %kSpace = gpuArray(kSpace);
    %Data.phase_mod_conj = gpuArray(Data.phase_mod_conj);
    %Data.sens_map_conj = gpuArray(Data.sens_map_conj);
    %new_img_z = gpuArray(new_img_z);
    %beta_sqrd = gpuArray(beta_sqrd);
%     gpuInfo = gpuDevice;
%     gpuSize = gpuInfo.AvailableMemory;
%     imSize  = numel(new_img_x);
    %     if imSize*para.Recon.no_comp > gpuSize*0.3
    %         para.Recon.type = [para.Recon.type,' less memory'];
    %     end
end

para.Cost = struct('fidelityNorm',[],'temporalNorm',[],'spatialNorm',[],'totalCost',[]);
fprintf(' Iteration       Cost       Step    Time(s) \n')
for iter_no = 1:para.Recon.noi
   
    t1 = tic;
    %fprintf('%.2f%%...',iter_no/noi_end*100);
%     fprintf(['Iteration = ' num2str(iter_no) '...']);

%%%%% fidelity term/temporal/spatial TV

    tic;
    [update_term, fidelity_norm] = compute_fidelity_yt_new(new_img_x,Data,para);
    para.CPUtime.fidelity(iter_no) = toc;
    
    tic;
    update_term = update_term + compute_3DtTV_yt(new_img_x,weight_tTV,beta_sqrd);
    para.CPUtime.tTV(iter_no) = toc;
    
    tic;
    update_term = update_term + compute_sTV_yt(new_img_x,weight_sTV,beta_sqrd);
    update_term = update_term + compute_sliceTV_yt(new_img_x,weight_slTV,beta_sqrd);
    para.CPUtime.sTV(iter_no) = toc;

%%%%% conjugate gradient
    tic;
    if iter_no > 1
        beta = update_term(:)'*update_term(:)/(update_term_old(:)'*update_term_old(:)+eps('single'));
        update_term = update_term + beta*update_term_old;
    end
    update_term_old = update_term; % clear update_term
    
%%%%% line search   

    %fidelity_update = compute_fidelity_for_line_search_yt(new_img_x,Data,para);
    
    para.Cost = Cost_STCR_3D(fidelity_norm, new_img_x, weight_sTV, weight_tTV, weight_slTV, para.Cost); clear fidelity_update
    step_size = line_search(new_img_x, update_term_old, Data, para);
    para.Recon.step_size(iter_no) = step_size;

    new_img_x = new_img_x + step_size * update_term_old;
    para.CPUtime.update(iter_no) = toc;

%%%%% plot&save part 

    if ifplot ==1
        showImage3D(new_img_x,para.Cost)
    end
    
    % break when step size too small or cost not changing too much
    if iter_no > 1 && para.Recon.break
        if step_size<1e-4 %|| abs(para.Cost.totalCost(end) - para.Cost.totalCost(end-1))/para.Cost.totalCost(end-1) < 1e-4
            break
        end
    end
    
    
fprintf(sprintf('%10.0f %10.2f %10.4f %10.2f \n',iter_no,para.Cost.totalCost(end),step_size,toc(t1)));
end

Image = squeeze(gather(new_img_x));
para = get_CPU_time(para);
fprintf(['Iterative STCR running time is ' num2str(para.CPUtime.interative_recon) 's' '\n'])
