function [kSpace,Kx,Ky,para] = PCA_3D_kSpace_yt(para)

tic;disp('Load k space data...');

PCA_dir          = para.dir.PCA_dir;
ifPCA            = para.Recon.ifPCA;

if ifPCA
    load(PCA_dir)
    if length(size(kSpace)) == 4
        para.Recon.Cart = 1;
    else
        para.Recon.Cart = 0;
    end
    disp('PCA result alreaty exists, skip PCA...');para.CPUtime.PCA = toc;toc;fprintf('\n');
else
    load(strcat(para.dir.load_kSpace_dir,para.dir.load_kSpace_name));
    if exist('PerfusionData','var')
        kSpace = PerfusionData; clear PerfusionData
    elseif exist('AIFData','var')
        kSpace = AIFData; clear AIFData
    end
    kSpace(isnan(kSpace)) = 0;
    if length(size(kSpace)) == 4
        para.Recon.Cart = 1;
        [sx,sy,sz,nc] = size(kSpace);
        nof = 1;
        Kx = 0;
        Ky = 0;
    else
        para.Recon.Cart = 0;
        [sx,sy,sz,nof,nc] = size(kSpace);
    end
    scale_kspace = 1e8;
    kSpace = kSpace*scale_kspace;
    para.CPUtime.load_kSpace = toc;toc;fprintf('\n');

%%%%% peforming PCA on coils

    disp('Peforming PCA on coils...');tic
    no_comp = para.NumberofPCAComponent;
    if no_comp == 'all'
        kSpace = single(kSpace);
        return
    end
    %data = double(kSpace);
    kSpace = reshape(kSpace,[sx*sy*sz*nof nc]);
    [coeff,~,~] = pca(kSpace); %principal component analysis
    kSpace = kSpace*coeff(:,1:no_comp);
    kSpace = reshape(kSpace,[sx sy sz nof no_comp]);
    clear coeff 
    save(PCA_dir,'kSpace','Kx','Ky');
    para.CPUtime.PCA = toc;toc;fprintf('\n');
end

if isfield(para,'FrameCutStart')
    kSpace = kSpace(:,:,:,para.FrameCutStart:para.FrameCutEnd,:);
    Kx = Kx(:,:,:,para.FrameCutStart:para.FrameCutEnd);
    Ky = Ky(:,:,:,para.FrameCutStart:para.FrameCutEnd);
    disp('Recon frames start:'),disp(para.FrameCutStart);
    disp('Recon frames end:'),disp(para.FrameCutEnd);
else
    disp('Recon all frames...');fprintf('\n');
end

