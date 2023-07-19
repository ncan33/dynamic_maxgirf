function [Image, para] = TFD_ADMM(Data, para)

disp('Performing iterative ADMM reconstruction...');
disp('Showing progress...')

disp_freq = 1;

ifplot         = para.setting.ifplot;
ifGPU          = para.setting.ifGPU;
weight_tTV     = para.Recon.weight_tTV;
weight_sTV     = para.Recon.weight_sTV;
beta_sqrd      = para.Recon.epsilon;
para.Recon.step_size = para.Recon.step_size(1);
% weight_l2      = para.Recon.weight_l2;


new_img_x = single(Data.first_est);

if isfield(Data,'sens_map')
    Data.sens_map_conj = conj(Data.sens_map);
end

if ifGPU
%    Data.kSpace        = gpuArray(Data.kSpace);
    new_img_x          = gpuArray(new_img_x);
    Data.sens_map      = gpuArray(Data.sens_map);
    Data.sens_map_conj = gpuArray(Data.sens_map_conj);
    if isfield(Data,'mask')
        Data.mask          = gpuArray(Data.mask);
    end
    if isfield(Data,'filter')
        Data.filter        = gpuArray(Data.filter);
    end
    beta_sqrd = gpuArray(beta_sqrd);
end

para.Cost = struct('fidelityNorm',[],'temporalNorm',[],'spatialNorm',[],'l2Norm',[],'totalCost',[]);
Cost = para.Cost;

fidelity = @(im) compute_fidelity_yt_new(im,Data,para);
temporal = @(im) compute_tTV_yt(im,weight_tTV,beta_sqrd);

V = zeros(size(new_img_x), class(new_img_x));
e = zeros(size(new_img_x), class(new_img_x));
weight_l2 = 10;

para.Recon.weight_l2 = weight_l2;
para.Cost.spatialNorm = 0;


fprintf(' Iteration       Cost       Step    Time(s) \n')
for iter_ADMM = 1:15
    Data.Ve = (V + e);
    dtVe = TFD(Data.Ve);
    dtVe = circshift(dtVe, [0, 0, 1]);
    
    for iter_CG = 1:15
        if mod(iter_CG,disp_freq) == 1 || iter_CG == 1 || disp_freq == 1
            t1 = tic;
        end
        
        %% fidelity term/temporal/spatial TV
        [update_term, fidelity_norm] = fidelity(new_img_x);

        dtm = TFD(new_img_x);
        
        dttm = TFD(dtm);
        dttm = circshift(dttm, [0, 0, 1]);
        
        tupdate = dttm - dtVe;
        
        update_term = update_term + weight_l2 * tupdate;
        
        %% conjugate gradient
        if iter_CG > 1
            beta = update_term(:)'*update_term(:)/(update_term_old(:)'*update_term_old(:)+eps('single'));
            update_term = update_term + beta*update_term_old;
        end
        update_term_old = update_term; clear update_term
        
        %% line search
        para.Cost = Cost_TFD_ADMM(fidelity_norm, gather(new_img_x), Data.Ve, weight_l2, para.Cost);
        Cost = Cost_STCR(fidelity_norm, new_img_x, weight_sTV, weight_tTV, Cost);
        %
        %     clear fidelity_update
        step_size = line_search_TFD_ADMM(new_img_x, update_term_old, Data, para);
        para.Recon.step_size(iter_CG) = step_size;
        %
        new_img_x = new_img_x + step_size * update_term_old;
        
        %% plot&save part
        if ifplot == 1
            showImage(new_img_x, Cost)
        end
        
        %% stoping creteria
        if para.Recon.break && iter_CG > 1
            if step_size<1e-4 %|| abs(para.Cost.totalCost(end) - para.Cost.totalCost(end-1))/para.Cost.totalCost(end-1) < 1e-4
                break
            end
        end
        
        if mod(iter_CG,disp_freq) == 0 || iter_CG == para.Recon.noi
            para.Recon.time(iter_CG) = toc(t1);
            fprintf(sprintf('%10.0f %10.2f %10.4f %10.2f \n',iter_CG,para.Cost.totalCost(end),step_size,para.Recon.time(iter_CG)));
        end
        
    end
    
    %% update V
    dtm = TFD(new_img_x);
    V = shrink_soft(dtm - e, weight_tTV/weight_l2);
    
    %% update e
    e = e - (dtm - V);
    
end

Image = squeeze(gather(new_img_x));
para.Recon.time_total = sum(para.Recon.time);
fprintf(['Iterative reconstruction running time is ' num2str(para.Recon.time_total) 's' '\n'])