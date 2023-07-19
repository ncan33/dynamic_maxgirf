function [kSpace_all,RayPosition,para] = load_kSpace_SMS_test(para)

tic;disp('Load k space data...');

ifasym           = para.asymmetry;
nSMS             = para.Recon.nSMS;
PCA_dir          = para.Recon.PCA_dir;


    load(strcat(para.dir.load_kSpace_dir,para.dir.load_kSpace_name));

    if(size(kSpace,4)==1)
        kSpace = permute(kSpace,[1 2 3 5 4]);
    end
    keyboard
    kSpace()
    scale_kspace = 10/mean(abs(kSpace(:)));
    kSpace = kSpace*scale_kspace;
    [sx,sy,nc,sz,ns] = size(kSpace);
    
    if ~exist('RayPosition','var')
        RayPosition = 1:sy;
        RayPosition = repmat(RayPosition',[1 sz]);
    end
    para.CPUtime.load_kSpace = toc;toc;fprintf('\n');
    
%%%%% peforming PCA on coils

    disp('Peforming PCA on coils...');tic
    no_comp = 8;
    data = double(permute(kSpace,[1 2 4 5 3]));
    data = reshape(data,[sx*sy*sz*ns nc]);
    [coeff,~,~] = pca(data); %principal component analysis
    compressed_data = data*coeff(:,1:no_comp);
    compressed_data = reshape(compressed_data,[sx sy sz ns no_comp]);
    compressed_data = permute(compressed_data,[1 2 3 5 4]);   
    kSpace = single(compressed_data);
    clear data compressed_data coeff score latent
    save(PCA_dir,'kSpace','RayPosition');
    para.CPUtime.PCA = toc;toc;fprintf('\n');


disp('Prepare k-space data, phase modulation and sampling angles...');tic

if isfield(para,'slice_pick')
    kSpace = kSpace(:,:,:,:,para.slice_pick);
end

if ifasym
    kSpace(para.AsymmetryRayStart:para.AsymmetryRayEnd,:,:,:,:) = 0;
end

RayPosition = reshape(RayPosition,[sy*sz 1]);
Ray_Index = find(RayPosition==1);

if nSMS == 3
    nor_SM = min(diff(Ray_Index)) - mod(min(diff(Ray_Index)),3); %This is very important. For 3 SMS
elseif length(Ray_Index) ~= 1
    nor_SM = min(diff(Ray_Index));
else
    nor_SM = 420;
end

if exist('FrameCutStart','var')
    Ray_Index = Ray_Index(FrameCutStart:FrameCutEnd);
    disp('Recon frames start:'),disp(FrameCutStart);
    disp('Recon frames end:'),disp(FrameCutEnd);
else
    disp('Recon all frames...');fprintf('\n');
end

nof = length(Ray_Index); %number of frames    

kSpace_all = reshape(kSpace,[sx,sy*sz,no_comp,1,ns]);

para.Recon.nor_SM = nor_SM;
para.Recon.nof = nof;
para.Recon.sy = sy;
para.Recon.sz = sz;




