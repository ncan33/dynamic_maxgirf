function Image = STCR_projections(Data,para)

disp('Performing conjugate gradient MBR...');
t1 = tic;

ifGPU          = para.setting.ifGPU;
ifplot         = para.setting.ifplot;
weight_sTV     = para.Recon.weight_sTV;
weight_tTV     = para.Recon.weight_tTV;
beta_sqrd      = para.Recon.epsilon;

%PD_abs = Data.PD_abs;

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
    Image = [];
    return
end

if isfield(Data,'phase_mod')
    Data.phase_mod_conj = conj(single(Data.phase_mod));
end
if isfield(Data,'sens_map')
    Data.sens_map_conj = conj(Data.sens_map);
end
%if isfield(Data,'first_guess')
%    new_img_x = Data.first_guess;
%else
%    new_img_x = single(Data.first_est);
%end
Image = Data.first_est;


if ifGPU
    %T1_temp = gpuArray(T1_temp);
    Image = gpuArray(Image);
    %Data.first_est = gpuArray(Data.first_est);
    Data.kSpace = gpuArray(Data.kSpace);
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

para.Cost = struct('fidelityNorm',[],'temporalNorm',[],'spatialNorm',[],'totalCost',[]);

for iter_no = 1:para.Recon.noi

    fidelity_update = compute_fidelity_yt_new(Image,Data,para);

    Image = reshape(Image,[288,288,30,10,3]);
    tTV_update = compute_tTV_yt(Image,weight_tTV,beta_sqrd);
    tTV_update = tTV_update + compute_3DtTV_yt(Image,weight_tTV/3,beta_sqrd);
    Image = reshape(Image,[288,288,300,1,3]);
    tTV_update = reshape(tTV_update,[288,288,300,1,3]);
    %tTV_update = compute_tTV_yt(Image,weight_tTV,beta_sqrd);
    sTV_update = compute_sTV_yt(Image,weight_sTV,beta_sqrd);
    update_term = fidelity_update + tTV_update + sTV_update;
    clear sTV_update tTV_update
    
    if iter_no > 1
        beta = update_term(:)'*update_term(:)/(update_term_old(:)'*update_term_old(:)+eps('single'));
        update_term = update_term + beta*update_term_old;
    end    
    update_term_old = update_term; clear update_term
    
    fidelity_update = compute_fidelity_for_line_search_yt(Image,Data,para);
    
    para.Cost = Cost_STCR(fidelity_update, Image, weight_sTV, weight_tTV, para.Cost); clear fidelity_update
    step_size = line_search(Image,update_term_old,Data,para);
    para.Recon.step_size(iter_no) = step_size;

    Image = Image + step_size * update_term_old;
    
    if ifplot ==1 %&& mod(iter_no,10) == 0
        showImage(Image,para.Cost)
    end


%{    
    if iter_no > 1 && para.Recon.break
        if step_size<1e-4 %|| abs(para.Cost.totalCost(end) - para.Cost.totalCost(end-1))/para.Cost.totalCost(end-1) < 1e-4
            break
        end
    end
%}


end

Image = squeeze(gather(Image));
para.PD_Cost = para.Cost;
para = rmfield(para,'Cost');
para.CPUtime.PD = toc(t1);

disp('PD images done.');toc(t1);


