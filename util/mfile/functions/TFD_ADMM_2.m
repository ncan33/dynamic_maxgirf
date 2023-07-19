function [Image, para] = TFD_ADMM_2(Data, para)

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

fidelity = @(im) compute_fidelity_yt_new(im,Data,para);
% temporal = @(im) compute_tTV_yt(im,weight_tTV,beta_sqrd);


V = zeros(size(new_img_x), class(new_img_x));
U = zeros([size(new_img_x), para.Recon.no_comp], class(new_img_x));
Z = zeros(size(new_img_x), class(new_img_x));

e1 = zeros(size(new_img_x), class(new_img_x));
e2 = zeros([size(new_img_x), para.Recon.no_comp], class(new_img_x));
e3 = zeros(size(new_img_x), class(new_img_x));

rho1 = 0.5;
rho2 = 0.5;
rho3 = 0.5;


%% Pre-calculate the inverse matrix
Sinv = rho2 + rho3;
Uinv = opt.U + opt.p2;
Tinv = invQT4d(opt);

para.Cost.spatialNorm = 0;


fprintf(' Iteration       Cost       Step    Time(s) \n')

for iter_no = 1:10
    if mod(iter_no,disp_freq) == 1 || iter_no == 1 || disp_freq == 1
        t1 = tic;
    end
    
    %% update Z
    temp = rho1 * TFD(V + e1) + rho3 * (new_img_x - e3);
    Z = ifftn(fftn(temp) ./ Tinv); 
    
    %% update V
    V = shrink(TFD(z) - e1, weight_tTV / rho1);
    
    %% update image
    new_img_x = (rho2 * sum(Data.sens_map_conj .* (U + e2), 4) + rho3 * (Z + e3)) / Sinv;
    
    %% update U
    temp = Data.sens_map .* new_img_x - e2;
    temp = rho2 * NUFFT.NUFFT(temp);
    temp = Data.kSpace + temp;
    
    U = NUFFT.NUFFT_adj(temp);
    
    %% update e1, e2, e3
    
    
    [update_term, fidelity_norm] = fidelity(new_img_x);
    
    %     dm = cat(3, new_img_x(:, :, 1), new_img_x, new_img_x(:, :, end));
    %     dm = diff(dm, 2, 3);
    
    dm = diff(new_img_x, 1, 3);
    dm(:, :, end+1) = new_img_x(:, :, 1) - new_img_x(:, :, end);
    
    ddm = diff(dm, 1, 3);
    ddm(:, :, end+1) = dm(:, :, 1) - dm(:, :, end);
    ddm = circshift(ddm, [0, 0, 1]);
    
    dVe = diff(Data.Ve, 1, 3);
    dVe(:, :, end+1) = Data.Ve(:, :, 1) - Data.Ve(:, :, end);
    %
    tupdate = ddm - dVe;
    
    update_term = update_term + weight_l2 * ddm;
    
    %     new_img_x = new_img_x + update_term * step_size;
    
    %% conjugate gradient
    if iter_no > 1
        beta = update_term(:)'*update_term(:)/(update_term_old(:)'*update_term_old(:)+eps('single'));
        update_term = update_term + beta*update_term_old;
    end
    update_term_old = update_term; clear update_term
    
    %% line search
    para.Cost = Cost_TFD_ADMM(fidelity_norm, gather(new_img_x), weight_tTV, Data.Ve, weight_l2, para.Cost);
    %
    %     clear fidelity_update
    step_size = line_search_TFD_ADMM(new_img_x,update_term_old,Data,para);
    para.Recon.step_size(iter_no) = step_size;
    %
    new_img_x = new_img_x + step_size * update_term_old;
    
    %% plot&save part
    if ifplot == 1
        showImage(new_img_x,para.Cost)
    end
    
    %% stoping creteria
    if para.Recon.break && iter_no > 1
        if step_size<1e-4 %|| abs(para.Cost.totalCost(end) - para.Cost.totalCost(end-1))/para.Cost.totalCost(end-1) < 1e-4
            break
        end
    end
    
    if mod(iter_no,disp_freq) == 0 || iter_no == para.Recon.noi
        para.Recon.time(iter_no) = toc(t1);
        fprintf(sprintf('%10.0f %10.2f %10.4f %10.2f \n',iter_no,para.Cost.totalCost(end),step_size,para.Recon.time(iter_no)));
    end
    
    
    
    dm = diff(new_img_x, 1, 3);
    dm(:, :, end + 1) = new_img_x(:, :, 1) - new_img_x(:, :, end);
    dm = dm - e;
    
    V = sign(dm) .* max(abs(dm) - weight_tTV / weight_l2, 0);
    
    e = e - (dm - V);
end

Image = squeeze(gather(new_img_x));
para.Recon.time_total = sum(para.Recon.time);
fprintf(['Iterative reconstruction running time is ' num2str(para.Recon.time_total) 's' '\n'])