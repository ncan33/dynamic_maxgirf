function Add_Dicoms_MSMS_SPGR_cine_perfusion(file,DicomFolder)


if ~isempty(strfind(DicomFolder,'/'))
    DicomFolderName = strfind(DicomFolder,'/');
    DicomFolderName = DicomFolder(DicomFolderName(end-1)+1:DicomFolderName(end)-1);
    Dicom_file_temp = dir([DicomFolder,'*.dcm']);
    Dicom_file_temp = [Dicom_file_temp(1).folder,'/',Dicom_file_temp(1).name];
else
    DicomFolderName = DicomFolder;
    Dicom_file_temp = DicomFolder;
end

load([file.folder,'/',file.name],'im_bin','para');
%im_bin = crop_half_FOV(abs(im_bin));
Image = permute(crop_half_FOV(im_bin),[1,2,4,3]);
Image = Image(:,:,1:3:end,:);
Image = abs(Image);

orintation = orintation_detection(sum(Image(:,:,1,:),4));
Image = orintate_image(Image,orintation);

[sx,sy,ns,nof] = size(Image);
% Image_sys = Image_sys(:,:,:,:);
% order = 1:nm:nm*ns;
% for i=ns:-1:2
%     order = [order;i:nm:nm*ns];
% end
% Image_sys = Image_sys(:,:,:,order);
% Nphase = size(para.Recon.bins,1);
% Ncycle = nof/Nphase;
% Image_sys = reshape(Image_sys,[sx,sy,Nphase,Ncycle,nm*ns]);


%im_bin = permute(im_bin,[1,2,3,5,4]);
%im_bin = reshape(im_bin,sx,sy,nof,nm*ns);
Image = Image/max(Image(:))*600;%scale

seriesNum = [file.name(4:8),'0'];
seriesDesc = para.dir.load_kSpace_name(1:end-14);
%seriesDesc = [seriesDesc,'F1000_T',num2str(para.weight_tTV*1000),'_S',num2str(para.weight_sTV*1000)];

parfor nSlice=1:ns
    DicomFolderLocalSys = ['ReconData/',DicomFolderName,'_cine_perfusion_slice_',num2str(nSlice,'%02.f')];
    mkdir(DicomFolderLocalSys)
    Image_one_slice_sys = squeeze(Image(:,:,nSlice,:));

%     Image_sys = squeeze(Image_one_slice(:,:,1,:));
%     Image_dia = squeeze(Image_one_slice(:,:,round((Nphase+1)/2),:));
    seriesNum_sys = [seriesNum,num2str(nSlice,'%02.f')];

    for nFrame = 1:nof
        outnamesys = [DicomFolderLocalSys,'/Img_',num2str(nFrame,'%03.f'),'.dcm'];
        addDicomHeader(Dicom_file_temp, Image_one_slice_sys(:,:,nFrame), seriesNum_sys, [seriesDesc,'cine_perfusion'], nFrame, outnamesys)
    end
end

end