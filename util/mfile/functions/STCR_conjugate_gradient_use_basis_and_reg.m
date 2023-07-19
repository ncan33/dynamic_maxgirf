function [Image,para] = STCR_conjugate_gradient_use_basis_and_reg(Data,para)
%[Image,para] = STCR_conjugate_gradient(Data,para)
disp('Performing iterative STCR reconstruction...');
disp('Showing progress...')

ifplot         = para.setting.ifplot;
ifGPU          = para.setting.ifGPU;
weight_tTV     = para.Recon.weight_tTV;
weight_sTV     = para.Recon.weight_sTV;
beta_sqrd      = para.Recon.epsilon;
para.Recon.step_size = para.Recon.step_size(1);

if isfield(para.Recon,'RF_frames')
    if sum(para.Recon.RF_frames)~=0
        RF_frames = para.Recon.RF_frames;
        PD_frames = para.Recon.PD_frames;
        Data.first_est = Data.first_est(:,:,RF_frames,:,:);
        Data.kSpace = Data.kSpace(:,:,RF_frames,:,:,:,:);
        if isfield(Data,'mask')
            Data.mask = Data.mask(:,:,RF_frames,:,:,:,:);
        end
        if isfield(Data,'phase_mod')
            Data.phase_mod = Data.phase_mod(:,:,RF_frames,:,:,:,:);
        end
        if isfield(Data,'N')
            if PD_frames(end) == 0
                PD_frames = 1:sum(PD_frames);
            end
            Data.N.S = Data.N.S(Data.N.sx_over.^2*PD_frames(end)+1:end,Data.N.siz(1)*Data.N.siz(2)*PD_frames(end)+1:end);
            Data.N.siz(3) = size(Data.first_est,3);
            %        Data.N.W = Data.N.W(:,:,RF_frames);
        end
    else
        Image = [];
        return
    end
elseif isfield(para.Recon,'PD_frames') && sum(para.Recon.PD_frames) == length(para.Recon.PD_frames)
    Image = [];
    return
end

if isfield(Data,'first_guess')
    new_img_x = Data.first_guess;   
else
    new_img_x = single(Data.first_est);
end

if isfield(Data,'phase_mod')
    Data.phase_mod_conj = conj(single(Data.phase_mod));
end

if isfield(Data,'sens_map')
    Data.sens_map_conj = conj(Data.sens_map);
end

if ifGPU
    Data.kSpace        = gpuArray(Data.kSpace);
    new_img_x          = gpuArray(new_img_x);
    Data.sens_map      = gpuArray(Data.sens_map);
    Data.sens_map_conj = gpuArray(Data.sens_map_conj);
    if isfield(Data,'mask')
        Data.mask          = gpuArray(Data.mask);
    end
    if isfield(Data,'filter')
        Data.filter        = gpuArray(Data.filter);
    end
%    Data.first_est = gpuArray(Data.first_est);
%    Data.phase_mod = gpuArray(Data.phase_mod);
%    Data.phase_mod_conj = gpuArray(Data.phase_mod_conj);
    beta_sqrd = gpuArray(beta_sqrd);
    if isfield(Data,'N')
        for i=1:length(Data.N)
            Data.N(i).S = gpuArray(Data.N(i).S);
            Data.N(i).Apodizer = gpuArray(Data.N(i).Apodizer);
            Data.N(i).W = gpuArray(Data.N(i).W);
        end
    end
    
    gpuInfo = gpuDevice;
    gpuSize = gpuInfo.AvailableMemory;
    imSize  = numel(new_img_x)*8;
%     if imSize*para.Recon.no_comp > gpuSize*0.3
        %para.Recon.type = [para.Recon.type,' less memory'];
%     end
end

para.Cost = struct('fidelityNorm',[],'temporalNorm',[],'spatialNorm',[],'totalCost',[]);
%new_img_x = new_img_x./Data.map;

fidelity = @(im) compute_fidelity_yt_new(im,Data,para);
spatial  = @(im) compute_sTV_yt(im,weight_sTV,beta_sqrd);
temporal = @(im) compute_tTV_yt(im,weight_tTV,beta_sqrd);




[sx,sy,nof,~,nSMS] = size(new_img_x);
nbasis = size(Data.basis{1},2);
for iter_no = 1:para.Recon.noi

    if mod(iter_no,10) == 1
        t1 = tic;
    end

    for i=1:nSMS
        Image_temp = new_img_x(:,:,:,:,i);
        Image_temp = reshape(Image_temp,[sx*sy,nof]);
        Image_temp = Data.basis{i}'*Image_temp.';
        
        Image_temp = Image_temp.';
        Image_temp = reshape(Image_temp,[sx,sy,nbasis]);
        
        [Model_Image] = T1_fitting(Image_temp,Data.Dic);
        Model_Image = reshape(Model_Image,[sx*sy,nbasis]);
        Model_Image = Model_Image.';
        
        Model_Image = Data.basis{i}*Model_Image;
        Model_Image = Model_Image.';
        Model_Image = reshape(Model_Image,[sx,sy,nof]);
        
        Image_temp = reshape(Image_temp,[sx*sy,nbasis]);
        Image_temp = Image_temp.';
        
        Image_temp = Data.basis{i}*Image_temp;
        Image_temp = Image_temp.';
        Image_temp = reshape(Image_temp,[sx,sy,nof]);
        
        
        new_img_x(:,:,:,:,i) = Image_temp;
        Model_Image_all(:,:,:,:,i) = Model_Image;
    end
    
keyboard
    
%     for i=1:nSMS
%         for j=1:nof
%             [~,x,y] = Reg_GS_tv(abs(new_img_x(:,:,j,1,i)),abs(Model_Image_all(:,:,j,1,i)),1,50);
%             new_img_x(:,:,j,1,i) = interp2(new_img_x(:,:,j,1,i),y,x);
%         end
%     end
            
%%%%% fidelity term/temporal/spatial TV

    tic; 
    [update_term,fidelity_norm] = fidelity(new_img_x);
    para.CPUtime.fidelity(iter_no) = toc;
    
    tic;
    update_term = update_term + temporal(new_img_x);
    para.CPUtime.tTV(iter_no) = toc;
    
    tic;
    update_term = update_term + spatial(new_img_x);
    para.CPUtime.sTV(iter_no) = toc;
    
    update_term = update_term + (Model_Image_all - new_img_x)*0.8;
    %update_term = update_term + (LowRank_yt(new_img_x) - new_img_x)*0.1;
    if isfield(para.Recon,'bins')
        update_term = update_term + compute_tTV_bins(new_img_x,weight_tTV,beta_sqrd,para.Recon.bins)*.5;
    end

%%%%% conjugate gradient
    tic;
    if iter_no > 1
        beta = update_term(:)'*update_term(:)/(update_term_old(:)'*update_term_old(:)+eps('single'));
        update_term = update_term + beta*update_term_old;
    end
    update_term_old = update_term; clear update_term
    
%%%%% line search    
    
    %fidelity_update = compute_fidelity_for_line_search_yt(new_img_x,Data,para);
    
    para.Cost = Cost_STCR(fidelity_norm, new_img_x, weight_sTV, weight_tTV, para.Cost); clear fidelity_update
    step_size = line_search(new_img_x,update_term_old,Data,para);
    para.Recon.step_size(iter_no) = step_size;

    new_img_x = new_img_x + step_size * update_term_old;
    para.CPUtime.update(iter_no) = toc;

%%%%% plot&save part 

    if ifplot == 1
        showImage(new_img_x,para.Cost)
    end

%%%%% break when step size too small or cost not changing too much

    if para.Recon.break && iter_no > 1
        if step_size<1e-4 %|| abs(para.Cost.totalCost(end) - para.Cost.totalCost(end-1))/para.Cost.totalCost(end-1) < 1e-4
            break
        end
    end
    
    if mod(iter_no,10) == 0
        fprintf(['Iteration = ' num2str(iter_no) '...']);
        toc(t1);
    end
end

Image = squeeze(gather(new_img_x));
para = get_CPU_time(para);
fprintf(['Iterative STCR running time is ' num2str(para.CPUtime.interative_recon) 's' '\n'])
end


function [Image,T1_map_all] = T1_fitting(Image,Dic)

nof = size(Image,3);

siz = [8*150/nof,nof/8];

fa = 0.4:0.05:1.2;
T1 = 100:1:2000;

if length(vec(Dic)) < 1000
    InversionTime = Dic;
    clear Dic
    for iT1=1:length(T1)
        for ifa = 1:length(fa)
            Dic(:,iT1,ifa) = Bloch_2_IR_1_SLG(InversionTime,0.92,T1(iT1),siz,flip_angle*fa(ifa),InversionTime(6)*1000-150);
        end
    end
    Dic = Dic(1:nof,:);
end
Dic_norm = sqrt(sum(Dic.^2));
%Dic_norm = min(Dic,[],1);
Dic = Dic./Dic_norm;
Dic = permute(Dic,[3,1,2]);
Dic = Dic(:,1:size(Image,3),:);

Image = squeeze(Image);
im_0 = abs(Image);
% non local mean filtering
% 
% im_0 = gather(im_0);
% for i=1:nof
%     for j=1:size(im_0,4)
%         im_0(:,:,i,j) = imnlmfilt(im_0(:,:,i,j),'DegreeOfSmoothing',0.1);
%     end
% end
% im_0 = gpuArray(im_0);
%%
im_p = angle(Image);
im_p = exp(1i.*im_p);
% im_p = im_p + compute_sTV_yt(im_p,0.1,eps('single'));


im_norm = sqrt(sum(im_0.^2,3));
%im_norm = min(im_0,[],3);
im_0 = im_0./im_norm;

siz = size(im_0);
if length(siz)<4
    siz(4) = 1;
end

%% fit
% T1_map_all = abs(zeros(siz(1),siz(2),siz(4),'like',Image));

tic
for i=1:siz(4)
    for j=1:siz(2)
        im = im_0(:,j,:,i);
        im = reshape(im,[siz(1),siz(3)]);
        
        d0 = im - abs(Dic);
        %d0(:,[1],:) = 0;
        %d0(:,1) = d0(:,1)/2;
        %d0(:,6) = d0(:,6)/2;
%        d0(:,[1:3,8:10,11:13,18:20,21:23,28:30,31:33,38:40,41:43,48:50],:,:) = [];
        d0 = sqrt(squeeze(sum(d0.^2,2)));
        
        [~,idx_T1] = min(min(d0,[],3),[],2);
        
        idx_map_all(:,j,i) = idx_T1;
    end
end  
toc

%

[T1_map_all,FA_map_all] = ind2sub([length(T1),length(fa)],idx_map_all);

% spatial constraint

% im_norm = im_norm + compute_sTV_yt(im_norm,3,1e-7);
% T1_map_all = T1_map_all + compute_sTV_yt(T1_map_all,20,1e-7);

%T1_map_all = round(T1_map_all);
%T1_map_all(T1_map_all>length(T1)) = length(T1);
%T1_map_all(T1_map_all<1) = 1;

idx_map_all = sub2ind([length(T1),length(fa)],T1_map_all,FA_map_all);
MBI = Dic(:,1:nof,round(idx_map_all));

T1_map_all = T1(T1_map_all);
FA_map_all = fa(FA_map_all);


MBI = reshape(MBI,[nof,siz([1,2,4])]);
MBI = permute(MBI,[2,3,1,4]);
MBI = abs(MBI).*im_norm;

%im_p = (mean(im_p,3)-im_p)*0.9+im_p;
%im_p = im_p + compute_tTV_yt(im_p,0.1,1e-7);
MBI = MBI(:,:,1:nof,:).*im_p;

%MBI = gpuArray(MBI);

%Image(73:216,73:216,:,:) = MBI; 
Image = permute(MBI,[1,2,3,5,4]);
end