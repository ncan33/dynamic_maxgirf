function [kSpace_all,RayPosition,para] = PCA_kSpace_unstreaking_yt(para)

tic;disp('Load k space data...');

%ifasym  = para.asymmetry;
nSMS    = para.Recon.nSMS;
PCA_dir = para.dir.PCA_dir;
%ifPCA   = para.Recon.ifPCA;
if ~strfind(PCA_dir,'unstreaking')
    PCA_dir = [PCA_dir(1:end-4),'_unstreaking.mat'];
end
ifPCA = ~isempty(dir(PCA_dir));

if ifPCA
    load(PCA_dir)
    [sx,sy,sz,no_comp,ns] = size(kSpace);
    disp('PCA result alreaty exists, skip PCA...');
    para.CPUtime.PCA = toc;toc;fprintf('\n');
else
    load(strcat(para.dir.load_kSpace_dir,para.dir.load_kSpace_name));
    if exist('AIFData','var')
        kSpace = AIFData; clear AIFData
        kSpace = permute(kSpace,[1 2 5 4 3]);
    end
    if(size(kSpace,4)==1)
        kSpace = permute(kSpace,[1 2 3 5 4]);
    end

    [sx,sy,nc,sz,ns] = size(kSpace);
    
    idx = sum(sum(sum(sum(kSpace,1),2),3),4)==0;
    %idx = idx==0;
    kSpace(:,:,:,:,idx) = [];
    ns = size(kSpace,5);
    
    %scale_kspace = 10^8;%sx^2/max(abs(kSpace(:)));
    %kSpace = kSpace*scale_kspace;
    
    if ~exist('RayPosition','var')
        RayPosition = 1:sy;
        RayPosition = repmat(RayPosition',[1 sz]);
    end
    kSpace_info.RayPosition = RayPosition;
    para.CPUtime.load_kSpace = toc;toc;fprintf('\n');
    
%%%%% peforming PCA on coils

    disp('Peforming PCA on coils...');tic
    if nc > para.NumberofPCAComponent
        kSpace = coil_unstreaking_and_compression_SMS(kSpace,kSpace_info);
        no_comp = para.NumberofPCAComponent;
        %data = double(permute(kSpace,[1 2 4 5 3]));
        %data = reshape(data,[sx*sy*sz*ns nc]);
        %[coeff,~,~] = pca(data); %principal component analysis
        %compressed_data = data*coeff(:,1:no_comp);
        %compressed_data = reshape(compressed_data,[sx sy sz ns no_comp]);
        %compressed_data = permute(compressed_data,[1 2 3 5 4]);
        %kSpace = single(compressed_data);
        %clear data compressed_data coeff score latent
        if exist('Kx','var')
            kSpace_info.kx = Kx;
            kSpace_info.ky = Ky;
        end        
        save(PCA_dir,'kSpace','kSpace_info');
    else
        no_comp = nc;
        kSpace = permute(kSpace,[1 2 4 3]);
    end

    para.CPUtime.PCA = toc;toc;fprintf('\n');
end

disp('Prepare k-space data, phase modulation and sampling angles...');tic

if isfield(para,'slice_pick') && ns~=1
    kSpace = kSpace(:,:,:,:,para.slice_pick);
    if size(kSpace_info.phase_mod,2)~=1
        kSpace_info.phase_mod = kSpace_info.phase_mod(:,para.slice_pick);
        para.phase_mod = para.phase_mod(:,para.slice_pick);
        para.angle_mod = para.angle_mod(:,para.slice_pick);
    end
    ns = length(para.slice_pick);
end

if isfield(para,'phase_correction')
    %kSpace = phase_correction_092817(kSpace);
    %kSpace = phase_correction_102017(kSpace);
    %kSpace = phase_correction_103117(kSpace);
    %kSpace = phase_correction_020218(kSpace);
    kSpace = phase_correction_031218(kSpace,kSpace_info.phase_mod);
end

%if ifasym
%    kSpace(para.AsymmetryRayStart:para.AsymmetryRayEnd,:,:,:,:) = 0;
%    kSpace_all = reshape(kSpace,[para.AsymmetryRayEnd,sy*sz,no_comp,1,ns]);
%else
kSpace_all = reshape(kSpace,[sx,sy*sz,no_comp,1,ns]);

%end
if exist('kSpace_info','var')
    RayPosition = kSpace_info.RayPosition;
end
RayPosition = reshape(RayPosition,[sy*sz 1]);
Ray_Index = find(RayPosition==1);
%RayPosition(Ray_Index(end):end) = [];
%Ray_Index(end) = [];

if nSMS == 3
    nor_SM = min(diff(Ray_Index)) - mod(min(diff(Ray_Index)),3); %This is very important. For 3 SMS
elseif length(Ray_Index) ~= 1
    nor_SM = min(diff(Ray_Index));
else
    nor_SM = 420;
end

%if isfield(para.Recon,'All_rays')
%    RayPosition(para.Recon.All_rays==0) = 0;
%end
if isfield(para.Recon,'SelectedFrames')
    Ray_Index = Ray_Index([para.Recon.SelectedFrames para.Recon.SelectedFrames(end)+1]);
    RayPosition(Ray_Index(end):end) = 0;
    RayPosition(1:Ray_Index(1)-1) = 0;
    %disp('Recon Selected Frames');
else
    %disp('Recon all frames...');fprintf('\n');
end
if isfield(para.Recon,'All_rays')
    all_rays = repmat(para.Recon.All_rays,[sy 1]);
    all_rays = all_rays(:);
    RayPosition(~all_rays) = 0;
end

%
nof = length(Ray_Index); %number of frames

[~,M] = max(abs(kSpace_all).^2,[],1);
M = mean(M(:));

if ~isfield(para,'kSpace_center')
    para.kSpace_center = M - 0.5;
elseif abs(para.kSpace_center - M + 0.5)>1
    %para.kSpace_center = M - 0.5;
end

para.Recon.nor_SM = nor_SM;
para.Recon.nof = nof;
para.Recon.sy = sy;
para.Recon.sz = sz;
para.Recon.sx = sx;
para.Recon.no_comp = no_comp;
para.Recon.image_size = [sx,sx];


