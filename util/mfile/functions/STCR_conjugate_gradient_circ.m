function [Image,para] = STCR_conjugate_gradient_circ(Data,para)
%--------------------------------------------------------------------------
%   [Image,para] = STCR_conjugate_gradient(Data,para)
%--------------------------------------------------------------------------
%   Solve MRI reconstruction problem using a conjugate gradient method.
%--------------------------------------------------------------------------
%   Inputs (for a 2D dynamic radial case):
%       - Data                      [structure] 
%           Data.kSpace             [sx, nor, nof, nc]
%           Data.sens_map           [1,  1,   1,   nc]
%           Data.first_est          [sx, sy,  nof]
%           Data.N                  [NUFFT structure]
%
%               'sx'    number of readout point along a ray
%               'sy'    for radial k-space, same as sx
%               'nor'   number of rays per time frame
%               'nof'   number of time frames
%               'nc'    number of coils
%           
%       - para                      [structure]
%           para.setting            [structure]
%               setting.ifplot      [0 or 1]
%               setting.ifGPU       [0 or 1]
%           para.Recon              [structure]
%               Recon.weight_tTV    [scalar]
%               Recon.weight_sTV    [scalar]
%               Recon.epsilon       [scalar]
%               Recon.step_size     [scalar]
%               Recon.noi           [positive integer]
%               Recon.type          [string]
%
%       - Data
%           Data.kSpace             measured k-space data "d"
%           Data.sens_map           sensitivity map
%           Data.first_est          initial estimation of "x": "A^H d"
%           Data.N                  NUFFT structure (see +NUFFT)
%
%       -para
%           para.setting.ifplot     display reconstruction process
%           para.setting.ifGPU      run function on a NVIDIA GPU
%           para.Recon.weight_tTV   "lambda_t"
%           para.Recon.weight_sTV   "lambda_s"
%           para.Recon.epsilon      "epsilon"
%           para.Recon.step_size    initial CG update step size
%           para.Recon.noi          number of iterations
%           para.Recon.type         reconstruction type see 
%                                   'compute_fidelity_ye_new'
%--------------------------------------------------------------------------
%   Output:
%       - Image     [sx, sy, nof, ...]
%       - para      [structure]
%
%       - Image     reconstructed images "m"
%--------------------------------------------------------------------------
%   A standard cost function it solves is the spatially and temporally
%   constrained reconstruction (STCR):
%
%   || Am - d ||_2^2 + lambda_t || TV_t m ||_1 + lambda_s || TV_s m ||_1
%
%   "A"         sampling matrix includes sensitivity maps, Fourier 
%               transform, and undersampling mask
%   "m"         image to be reconstructed
%   "d"         measured k-space data
%   ||.||_2^2   l2 norm
%   ||.||_1     l1 norm
%   "lambda_t"  temporal constraint weight
%   "lambda_s"  sparial constraint weight
%   TV_t        temporal total variation (TV) operator (finite difference)
%               sqrt( abs(m_t+1 - m_t)^2 + epsilon )
%   "epsilon"   small term to aviod singularity
%   TV_s        spatial TV operator
%               sqrt( abs(m_x+1 - m_x)^2 + abs(m_y+1 - m_y) + epsilon )
%--------------------------------------------------------------------------
%   Reference:
%       [1]     Acquisition and reconstruction of undersampled radial data 
%               for myocardial perfusion MRI. JMRI, 2009, 29(2):466-473.
%--------------------------------------------------------------------------
%   Author:
%       Ye Tian
%       E-mail: phye1988@gmail.com
%--------------------------------------------------------------------------

%% load reconstruction parameters and settings

fprintf([repmat('-', [1, 75]), '\n'])
disp('begin iterative STCR conjugate gradient reconstruction...');
fprintf([repmat('-', [1, 75]), '\n'])

disp_freq               = 1;
ifplot                  = para.setting.ifplot;
ifGPU                   = para.setting.ifGPU;
weight_tTV              = para.Recon.weight_tTV;
weight_sTV              = para.Recon.weight_sTV;
beta_sqrd               = para.Recon.epsilon;
para.Recon.step_size    = para.Recon.step_size(1);

%{
%lambda_vbm4d   = para.vbm4d;
%lambda_new    = para.recon.low_rank;
%bloc_x = 32;
%bloc_y = 32;
%tau = (0.01/5)*bloc_x;
%lambda_new = para.Recon.low_rank;
% VBM4D params

% intensity range of the video
%profile = 'np';      % V-BM4D parameter profile
%  'lc' --> low complexity
%  'np' --> normal profile
%do_wiener = 1;       % Wiener filtering
%   1 --> enable Wiener filtering
%   0 --> disable Wiener filtering
%sharpen = 1;         % Sharpening
%   1 --> disable sharpening
%  >1 --> enable sharpening
%deflicker = 1;       % Deflickering
%   1 --> disable deflickering
%  <1 --> enable deflickering
%verbose = 0;         % Verbose mode
%p_sigma     = 1;      % Noise standard deviation. it should be in the same
%}
if isfield(para.Recon,'RF_frames')
    if sum(para.Recon.RF_frames)~=0
        RF_frames = para.Recon.RF_frames;
        PD_frames = para.Recon.PD_frames;
        nof = size(Data.first_est,3);
        Data.first_est = Data.first_est(:,:,RF_frames(1:nof),:,:);
        Data.kSpace = Data.kSpace(:,:,RF_frames(1:nof),:,:,:,:);
        if isfield(Data,'mask')
            Data.mask = Data.mask(:,:,RF_frames(1:nof),:,:,:,:);
        end
        if isfield(Data,'phase_mod')
            Data.phase_mod = Data.phase_mod(:,:,RF_frames(1:nof),:,:,:,:);
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
    if ~isfield(Data, 'sens_map_conj')
        Data.sens_map_conj = conj(Data.sens_map);
    end
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
    if isfield(Data, 'off_res')
        Data.off_res.f_im = gpuArray(Data.off_res.f_im);
        Data.off_res.f_k = gpuArray(Data.off_res.f_k);
    end
    if isfield(Data, 'phase_mod')
        Data.phase_mod = gpuArray(Data.phase_mod);
        Data.phase_mod_conj = gpuArray(Data.phase_mod_conj);
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
    
%     gpuInfo = gpuDevice;
%     gpuSize = gpuInfo.AvailableMemory;
%     imSize  = numel(new_img_x)*8;
%     if imSize*para.Recon.no_comp > gpuSize*0.3
        %para.Recon.type = [para.Recon.type,' less memory'];
%     end
end

% initialize Cost
para.Cost = struct('fidelityNorm',[],'temporalNorm',[],'spatialNorm',[],'totalCost',[]);

% initialize function handels 
fidelity = @(im) compute_fidelity_yt_new(im,Data,para);
spatial  = @(im) compute_sTV_yt(im,weight_sTV,beta_sqrd);
temporal = @(im) compute_tTV_circ(im,weight_tTV,beta_sqrd);

fprintf(' Iteration       Cost       Step    Time(s) \n')
for iter_no = 1:para.Recon.noi

    if mod(iter_no,disp_freq) == 1 || iter_no == 1 || disp_freq == 1
        t1 = tic;
    end

%% update term (gradient of fidelity, temporal/spatial TV)

    tic; 
    [update_term,fidelity_norm] = fidelity(new_img_x);
    para.CPUtime.fidelity(iter_no) = toc;
    
    tic;
    update_term = update_term + temporal(new_img_x);
    para.CPUtime.tTV(iter_no) = toc;
    
    tic;
    update_term = update_term + spatial(new_img_x);
    para.CPUtime.sTV(iter_no) = toc;

%% conjugate gradient
    tic;
    if iter_no > 1
        beta = update_term(:)'*update_term(:)/(update_term_old(:)'*update_term_old(:)+eps('single'));
        update_term = update_term + beta*update_term_old;
    end
    update_term_old = update_term; clear update_term
    
%% line search    
    
    %fidelity_update = compute_fidelity_for_line_search_yt(new_img_x,Data,para);
    
    para.Cost = Cost_STCR(fidelity_norm, new_img_x, weight_sTV, weight_tTV, para.Cost); clear fidelity_update
    step_size = line_search(new_img_x,update_term_old,Data,para);
    para.Recon.step_size(iter_no) = step_size;

    new_img_x = new_img_x + step_size * update_term_old;
    para.CPUtime.update(iter_no) = toc;
%{
%% VBM4D

    tic;
    if(lambda_vbm4d~=0)
        vbm4d_term_update = VBM4D_ye(new_img_x);
        new_img_x = new_img_x + lambda_vbm4d*(vbm4d_term_update-new_img_x);
    end
    para.CPUtime.VBM4D(iter_no) = toc;
 
%% Low Rank

    if lambda_new~=0
        for i=1:size(new_img_x,5)
            lr_update(:,:,:,:,i) = low_rank_yt(new_img_x(:,:,:,:,i),bloc_x,bloc_y,tau);
        end

        new_img_x = new_img_x + lambda_new * (lr_update - new_img_x);
        clear lr_update
    end
%}
%% plot&save part 

    if ifplot == 1
        showImage(new_img_x,para.Cost)
    end

%% break when step size too small or cost not changing too much

    if para.Recon.break && iter_no > 1
        if step_size < 1e-5 %|| abs(para.Cost.totalCost(end) - para.Cost.totalCost(end-1))/para.Cost.totalCost(end-1) < 1e-4
            break
        end
    end
    
    if mod(iter_no,disp_freq) == 0 || iter_no == para.Recon.noi
        fprintf(sprintf('%10.0f %10.2f %10.4f %10.2f \n',iter_no,para.Cost.totalCost(end),step_size,toc(t1)));
    end
end

Image = gather(new_img_x);
para = get_CPU_time(para);
fprintf(['Iterative STCR running time is ' num2str(para.CPUtime.interative_recon) 's' '\n'])
fprintf([repmat('-', [1, 75]), '\n'])