function para = load_default_para(option)

para.setting.PCA = 1;% do PCA or not
para.NumberofPCAComponent = 8;

if strfind(option,'SMS')
    para.phase_mod = 1;
else
    para.phase_mod = 0;
end
if strfind(option,'golden')
    para.angle_mod = 1;
end
if strfind(option,'cSMS')
    para.cSMS = 1;
else
    para.cSMS = 0;
end
if contains(option,'GROG')
    para.Recon.interp_method = 'GROG';
    para.over_sampling = 1;
    para.core_size = [1,1];
end
if contains(option,'NUFFT')
    para.Recon.interp_method = 'NUFFT';
    para.over_sampling = 1.5;
    para.core_size = [6,6];
end
if contains(option,'GPU')
    para.setting.ifGPU = 1;
else
    para.setting.ifGPU = 0;
end

para.setting.save_frequency = 20000;
para.setting.debug = 1;
para.setting.ifplot = 1;

%%% ifNUFFT = 0, then a method for interpolation has to be chosen. 
%%% 'nn' for nearest neighbor, 'grid3' and 'GROG'
para.ifNUFFT = 0;
%para.interp_method = 'GROG';
%para.over_sampling_factor = 1;

para.setting.BacktrackingLineSearch = 1;
para.setting.FirstEstimationA3J = 0;

para.Recon.epsilon = eps('single');

para.asymmetry = 0;

if strfind(option,'360')
    para.dataAngle = 360;
else
    para.dataAngle = 180;
end
%para.ifGPU = 0;

para.vbm4d = 0;

para.weight_sTV = 0.001;
para.weight_tTV = 0.02;

if strfind(option,'3D')
    para.weight_sTV = 0.0005;
    para.weight_sliceTV = 0.0005;
end

para.Recon.step_size = 2; %0.5 is good for NUFFT
para.Recon.crop_half_FOV = 1;
para.image_orintation = 1;
para.noi = 150;
%para.kSpace_center = 145.2;

para.dir.save_recon_img_mat_dir = strcat(pwd,'/ReconData/mat_files/');