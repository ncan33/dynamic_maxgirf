function para = prepare_para(para)
fprintf('Loading parameters...');tic

matObj = matfile([para.dir.load_kSpace_dir,para.dir.load_kSpace_name]);
varlist = who(matObj,'kSpace_info');
if ~isempty(varlist)
    load([para.dir.load_kSpace_dir,para.dir.load_kSpace_name],'kSpace_info')
    if isfield(kSpace_info,'phase_mod')
        para.phase_mod = kSpace_info.phase_mod;
    end
    if isfield(kSpace_info,'angle_mod')
        para.angle_mod = kSpace_info.angle_mod;
    end
    if isfield(kSpace_info,'ray_mod')
        ray_mod = kSpace_info.ray_mod(:,1);
        ray_mod = reshape(ray_mod,kSpace_info.Nr,kSpace_info.Nf);
        ray_mod = ray_mod(1,:);
        para.Recon.PD_frames = ray_mod == 0;
        if sum(para.Recon.PD_frames) ~= length(ray_mod) && ~isfield(para.Recon,'RF_frames')
            para.Recon.RF_frames = ray_mod == ray_mod(sum(para.Recon.PD_frames)+1);
            para.Recon.All_rays = para.Recon.PD_frames + para.Recon.RF_frames;
        end
    end
    if isfield(kSpace_info,'line_type')
        ray_mod = kSpace_info.line_type;
        para.Recon.PD_frames = false(1,max(kSpace_info.frames));
        para.Recon.PD_frames(1:kSpace_info.frames(find(diff(ray_mod==0)))) = true;
        para.Recon.RF_frames = ~para.Recon.PD_frames;
    end
end

para.time = datestr(clock,'yymmdd_hhMMSS');

if isfield(para.dir,'continue_dir')
    para.Recon.ifContinue = 1;
else
    para.Recon.ifContinue = 0;
end

%if para.Recon.ifContinue
%    para.Recon.noi_start = para.noi_old;
%    para.Recon.noi_end = para.noi_old + para.noi_add;
%else
%    para.Recon.noi_start = 0;
%    para.Recon.noi_end = para.noi;
%end

switch para.Recon.interp_method
    case 'NUFFT'
    para.ifNUFFT = 1;
end
%ifNUFFT = para.ifNUFFT;
%ifA3J = para.FirstEstimationA3J;

if para.phase_mod == 0
    nSMS = 1;
else
    nSMS = 3;
end
para.Recon.nSMS = nSMS;
MID_i = strfind(para.dir.load_kSpace_name,'MID');
para.MID = para.dir.load_kSpace_name(MID_i+5:MID_i+7);
MID = para.MID;

disp('RawData:')
disp([para.dir.load_kSpace_dir para.dir.load_kSpace_name])

try
    load([para.dir.load_kSpace_dir,'*',MID,'*data_info*'])
    para.data_info = data_info;
catch
    disp('No data_info')
end
%selected_rays_start    = para.selected_rays_start;
%selected_rays_end      = para.selected_rays_end;



%weight_tTV          = para.weight_tTV;
%weight_sTV          = para.weight_sTV;

%ifGPU = para.ifGPU;
%{
if isempty(dir([pwd,'/SensitivityMap/']))
    mkdir([pwd,'/SensitivityMap/']);
end
%}
%if ~isfield(para.Recon,'sens_map_dir')
%    if ifNUFFT==1
%        sens_map_dir = strcat('/SensitivityMap/s_map_NUFFT_MID_',int2str(MID),'.mat');
%    else
%        interp_method = para.interp_method;
%        sens_map_dir = strcat('/SensitivityMap/s_map_PreInterp_',interp_method, '_MID_',int2str(MID),'.mat');
%    end
%    para.Recon.sens_map_dir = sens_map_dir;
%else
%    sens_map_dir = para.Recon.sens_map_dir;
%end



kSpace_name = char(para.dir.load_kSpace_name);

name = '';
if ~isempty(strfind(kSpace_name,'perfusion')) || ~isempty(strfind(kSpace_name,'Perfusion')) || ~isempty(strfind(kSpace_name,'PERFUSION'))
    name = strcat(name,'Perfusion_');
end
if ~isempty(strfind(kSpace_name,'Rest')) || ~isempty(strfind(kSpace_name,'rest')) || ~isempty(strfind(kSpace_name,'REST'))
    name = strcat(name,'Rest_');
elseif ~isempty(strfind(kSpace_name,'Stress')) || ~isempty(strfind(kSpace_name,'stress')) || ~isempty(strfind(kSpace_name,'STRESS'))
    name = strcat(name,'Stress_');
end
if ~isempty(strfind(kSpace_name,'AIF'))
    name = strcat(name,'AIF_');
end
if ~isempty(strfind(kSpace_name,'3D'))
    name = strcat(name,'3D_');
end
if ~isempty(strfind(kSpace_name,'echo_1'))
    name = strcat(name,'echo_1_');
elseif ~isempty(strfind(kSpace_name,'echo_2'))
    name = strcat(name,'echo_2_');
elseif ~isempty(strfind(kSpace_name,'echo_3'))
    name = strcat(name,'echo_3_');
end
if ~isempty(strfind(kSpace_name,'diastole')) || ~isempty(strfind(kSpace_name,'Diastole'))
    name = strcat(name,'diastole_');
elseif ~isempty(strfind(kSpace_name,'systole')) || ~isempty(strfind(kSpace_name,'Systole'))
    name = strcat(name,'systole_');
end

if isfield(para,'name_add')
    name = strcat(name,para.name_add,'_');
    para = rmfield(para,'name_add');
end

%name = [num2str(para.over_sampling),'_',num2str(para.core_size(1)),'x',num2str(para.core_size(2)),'_',name];

para.dir.save_recon_img_name= strcat('SET',num2str(para.set,'%05.f_'),'MID',num2str(MID),'_',para.Recon.interp_method,'_',name);

para.dir.save_recon_img_name = strcat(para.dir.save_recon_img_name,para.time,'.mat');
%{
if ~isfield(para.Recon,'sens_map_dir')
    sens_map_dir = strcat('/SensitivityMap/s_map_MID',num2str(MID),'_',name(1:end-1),'.mat');
    para.Recon.sens_map_dir = sens_map_dir;
else
    sens_map_dir = para.Recon.sens_map_dir;
end

sens_map_precomputed = dir([pwd,sens_map_dir]);
noSensMap = isempty(sens_map_precomputed);
para.Recon.noSensMap = noSensMap;
%}
%para.Recon.save_dir = strcat(para.Recon.save_dir,num2str(weight_tTV,'_T%#.4f'),num2str(weight_sTV,'_S%#.4f_'));


%save_recon_img_name = strcat('rays_',num2str(selected_rays_start),'_',num2str(selected_rays_end));
%if ifA3J ==1
%    save_recon_img_name = strcat('A3J_',save_recon_img_name);
%end

%if ifNUFFT==1
%    save_recon_img_name = strcat('NUFFT_',save_recon_img_name);
%else
%    interp_method = para.interp_method;
%    switch interp_method
%        case 'grid3'
%            save_recon_img_name = strcat('grid3_OIF_',num2str(para.over_sampling_factor),'_',save_recon_img_name);
%        case 'nn'
%            save_recon_img_name = strcat('nn_',save_recon_img_name);
%        case 'GROG'
%            save_recon_img_name = strcat('GROG_',save_recon_img_name);
%    end
%end
%if ifGPU ==1
%    save_recon_img_name = strcat('GPU_',save_recon_img_name);
%end

if isempty(dir(para.dir.save_recon_img_mat_dir))
    mkdir(para.dir.save_recon_img_mat_dir);
end

if isempty(dir([pwd,'/RawData/']))
    mkdir([pwd,'/RawData/']);
end

kSpace_data_dir  = para.dir.load_kSpace_dir;
kSpace_data_name = para.dir.load_kSpace_name;
if isfield(para.setting,'unstreaking')
    if para.setting.unstreaking
        PCA_name = strcat(kSpace_data_name(1:end-4),'_PCA_unstreaking.mat');
    else
        PCA_name = strcat(kSpace_data_name(1:end-4),'_PCA.mat');
    end
else
    PCA_name = strcat(kSpace_data_name(1:end-4),'_PCA.mat');
end
PCA_dir_oringinal = [kSpace_data_dir,PCA_name];
if isempty(dir(PCA_dir_oringinal))
    para.dir.PCA_dir = [pwd,'/RawData/',PCA_name];
else
    para.dir.PCA_dir = PCA_dir_oringinal;
end

PCA_precomputed = dir(para.dir.PCA_dir);
para.Recon.ifPCA = ~isempty(PCA_precomputed);

matObj = matfile(para.dir.PCA_dir);
varlist = who(matObj,'kSpace_info');
if ~isempty(varlist)
    load(para.dir.PCA_dir,'kSpace_info')
    if isfield(kSpace_info,'phase_mod')
        para.phase_mod = kSpace_info.phase_mod;
    end
    if isfield(kSpace_info,'angle_mod')
        para.angle_mod = kSpace_info.angle_mod;
    end
    if isfield(kSpace_info,'ray_mod')
        ray_mod = kSpace_info.ray_mod(:,1);
        ray_mod = reshape(ray_mod,kSpace_info.Nr,kSpace_info.Nf);
        ray_mod = ray_mod(1,:);
        para.Recon.PD_frames = ray_mod == 0;
        if sum(para.Recon.PD_frames) ~= length(ray_mod)
            para.Recon.RF_frames = ray_mod == ray_mod(sum(para.Recon.PD_frames)+1);
            para.Recon.All_rays = para.Recon.PD_frames + para.Recon.RF_frames;
        end
    end
    if isfield(kSpace_info,'line_type')
        ray_mod = kSpace_info.line_type;
        if isfield(kSpace_info,'frames')
            para.Recon.PD_frames = false(1,max(kSpace_info.frames));
            para.Recon.PD_frames(1:kSpace_info.frames(find(diff(ray_mod==0)))) = true;
        else
            para.Recon.PD_frames = false(1,length(kSpace_info.MeasurementTime_seconds));
            para.Recon.PD_frames(1:find(diff(ray_mod))/kSpace_info.RadialViews) = true;
        end
        para.Recon.RF_frames = ~para.Recon.PD_frames;
    end
end

if ~isfield(para.Recon,'PD_frames') && exist('kSpace_info')
    if isfield(kSpace_info,'MeasurementTime_seconds')
        para.Recon.PD_frames = false(1,length(kSpace_info.MeasurementTime_seconds));
        para.Recon.PD_frames(1:kSpace_info.ProtonDensityScans) = true;
        para.Recon.RF_frames = ~para.Recon.PD_frames;
    end
end

para.dir.cart_kSpace_dir = strcat(kSpace_data_dir,kSpace_data_name(1:end-4),'_cart.mat');
interp_precomputed = dir(para.dir.cart_kSpace_dir);
para.Recon.ifInterp = ~isempty(interp_precomputed);
if exist('kSpace_info')
    para.kSpace_info = kSpace_info;
end

para.dir.save_recon_img_mat_dir = [para.dir.save_recon_img_mat_dir,'MID',num2str(MID),'/'];
if isempty(dir(para.dir.save_recon_img_mat_dir))
    mkdir(para.dir.save_recon_img_mat_dir)
end

para.CPUtime.load_para_time = toc;toc;fprintf('\n');
end
