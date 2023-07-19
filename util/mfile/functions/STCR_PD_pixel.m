function [Image_PD,para] = STCR_PD_pixel(Data,para)

disp('Performing conjugate gradient STCR on PD...');
t1 = tic;

ifGPU          = para.setting.ifGPU;
ifplot         = para.setting.ifplot;
weight_sTV     = para.Recon.weight_sTV;
weight_tTV     = para.Recon.weight_tTV;
beta_sqrd      = para.Recon.epsilon;
%para.Recon.step_size = para.Recon.step_size(1);
%step_size_0    = para.Recon.step_size;

if isfield(para.Recon,'PD_frames') && sum(para.Recon.PD_frames)~=0
    PD_frames = para.Recon.PD_frames;
    if PD_frames(end) == 0
        PD_frames = 1:sum(PD_frames);
    end
    Data.first_est = Data.first_est(:,:,PD_frames,:,:);
    Data.kSpace = Data.kSpace(:,:,PD_frames,:,:,:,:);
    if isfield(Data,'mask')
        Data.mask = Data.mask(:,:,PD_frames,:,:,:,:);
    end
    if isfield(Data,'phase_mod')
        Data.phase_mod = Data.phase_mod(:,:,PD_frames,:,:,:,:);
    end
    if isfield(Data,'N')
        Data.N.S = Data.N.S(1:Data.N.sx_over.^2*PD_frames(end),1:Data.N.siz(1)*Data.N.siz(2)*PD_frames(end));
        Data.N.siz(3) = size(Data.first_est,3);
        Data.N.W = Data.N.W(:,:,PD_frames);
    end
else
    Image_PD = [];
    return
end


if isfield(Data,'phase_mod')
    Data.phase_mod_conj = conj(single(Data.phase_mod));
end
if isfield(Data,'sens_map')
    Data.sens_map_conj = conj(Data.sens_map);
end

if isfield(Data,'first_guess')
    new_img_x = Data.first_guess;
    noi_start = length(para.Recon.step_size)+1;
    para.Recon.noi = para.Recon.noi + noi_start;
    para.Cost = para.PD_Cost;
else
    new_img_x = Data.first_est;
    noi_start = 1;
    para.Cost = struct('fidelityNorm',[],'temporalNorm',[],'spatialNorm',[],'totalCost',[]);
end


if ifGPU
    new_img_x = gpuArray(new_img_x);
    Data.first_est = gpuArray(Data.first_est);
    Data.sens_map = gpuArray(Data.sens_map);
    Data.sens_map_conj = gpuArray(Data.sens_map_conj);
    beta_sqrd = gpuArray(beta_sqrd);
    if isfield(Data,'N')
        for i=1:length(Data.N)
            Data.N(i).S = gpuArray(Data.N(i).S);
            Data.N(i).Apodizer = gpuArray(Data.N(i).Apodizer);
            Data.N(i).W = gpuArray(Data.N(i).W);
        end
    end
end



for iter_no = noi_start:para.Recon.noi

    fidelity_update = compute_fidelity_yt_new(new_img_x,Data,para);
    sTV_update = compute_sTV_yt(new_img_x,weight_sTV,beta_sqrd);
    tTV_update = compute_tTV_pixel(new_img_x,weight_tTV,beta_sqrd,Data.Motion);
    update_term = fidelity_update + tTV_update + sTV_update;
    clear sTV_update tTV_update
    
    if iter_no > noi_start
        beta = update_term(:)'*update_term(:)/(update_term_old(:)'*update_term_old(:)+eps('single'));
        update_term = update_term + beta*update_term_old;
    end
    
    update_term_old = update_term; clear update_term
    
    fidelity_update = compute_fidelity_for_line_search_yt(new_img_x,Data,para);
    
    para.Cost = Cost_STCR_pixel(fidelity_update, new_img_x, Data.Motion, weight_sTV, weight_tTV, para.Cost); clear fidelity_update
    step_size = line_search_pixel(new_img_x,update_term_old,Data,para);
    para.Recon.step_size(iter_no) = step_size;

    new_img_x   = new_img_x + step_size * update_term_old;
    
    if ifplot ==1 %&& mod(iter_no,10) == 0
        showImage(new_img_x,para.Cost)
    end

    if iter_no > 1 && para.Recon.break
        if step_size<1e-4 %|| abs(para.Cost.totalCost(end) - para.Cost.totalCost(end-1))/para.Cost.totalCost(end-1) < 1e-4
            break
        end
    end

end

%showImage(new_img_x,para.Cost)
Image_PD = squeeze(gather(abs(new_img_x)));
%Image_PD = mean(abs(Image_PD),3);
%Data.first_est(:,:,PD_frames,:,:) = [];
%Data.mask(:,:,PD_frames,:,:,:,:) = [];
para.PD_Cost = para.Cost;
para.PD_step_size = para.Recon.step_size;
%para.Recon.step_size = step_size_0;
para = rmfield(para,'Cost');
para.CPUtime.PD = toc(t1);

disp('PD images done.');toc(t1);