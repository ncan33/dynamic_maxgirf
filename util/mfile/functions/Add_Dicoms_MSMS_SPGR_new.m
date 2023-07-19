function Add_Dicoms_MSMS_SPGR_new(file,DicomFolder)


if ~isempty(strfind(DicomFolder,'/'))
    DicomFolderName = strfind(DicomFolder,'/');
    DicomFolderName = DicomFolder(DicomFolderName(end-1)+1:DicomFolderName(end)-1);
    Dicom_file_temp = dir([DicomFolder,'*.dcm']);
    Dicom_file_temp = [Dicom_file_temp(1).folder,'/',Dicom_file_temp(1).name];
else
    DicomFolderName = DicomFolder;
    Dicom_file_temp = DicomFolder;
end

load([file.folder,'/',file.name],'Image_sys','Image_dia','para')
%im_bin = crop_half_FOV(abs(im_bin));
if size(Image_sys,4)<size(Image_sys,3)
    Image_sys = permute(Image_sys,[1,2,4,3]);
end
if size(Image_dia,4)<size(Image_dia,3)
    Image_dia = permute(Image_dia,[1,2,4,3]);
end
if size(Image_sys,3)>12
    Image_sys = Image_sys(:,:,1:3:end,:);
    Image_dia = Image_dia(:,:,1:3:end,:);
end
orintation = orintation_detection(sum(Image_sys(:,:,1,:),4));

Image_sys = orintate_image(Image_sys,orintation);
Image_dia = orintate_image(Image_dia,orintation);

[sx,sy,ns,nof] = size(Image_sys);
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
Image_sys = Image_sys/max(Image_sys(:))*600;%scale

seriesNum = file.name(4:8);
seriesDesc = para.dir.load_kSpace_name(1:end-14);
%seriesDesc = [seriesDesc,'F1000_T',num2str(para.weight_tTV*1000),'_S',num2str(para.weight_sTV*1000)];

parfor nSlice=1:ns
    DicomFolderLocalSys = ['ReconData/',DicomFolderName,'_sys_slice_',num2str(nSlice,'%02.f')];
    mkdir(DicomFolderLocalSys)
    DicomFolderLocalDia = ['ReconData/',DicomFolderName,'_dia_slice_',num2str(nSlice,'%02.f')];
    mkdir(DicomFolderLocalDia)
    Image_one_slice_sys = squeeze(Image_sys(:,:,nSlice,:));
    Image_one_slice_dia = squeeze(Image_dia(:,:,nSlice,:));
%     Image_sys = squeeze(Image_one_slice(:,:,1,:));
%     Image_dia = squeeze(Image_one_slice(:,:,round((Nphase+1)/2),:));
    seriesNum_sys = [seriesNum,num2str(nSlice,'%02.f')];
    seriesNum_dia = [seriesNum,num2str(nSlice+50,'%02.f')];
    for nFrame = 1:nof
        outnamesys = [DicomFolderLocalSys,'/Img_',num2str(nFrame,'%03.f'),'.dcm'];
        outnamedia = [DicomFolderLocalDia,'/Img_',num2str(nFrame,'%03.f'),'.dcm'];
        addDicomHeader(Dicom_file_temp, Image_one_slice_sys(:,:,nFrame), seriesNum_sys, [seriesDesc,'sys'], nFrame, outnamesys)
        addDicomHeader(Dicom_file_temp, Image_one_slice_dia(:,:,nFrame), seriesNum_dia, [seriesDesc,'dia'], nFrame, outnamedia)
    end
end

end