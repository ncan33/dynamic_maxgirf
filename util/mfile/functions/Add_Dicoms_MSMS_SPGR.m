function Add_Dicoms_MSMS_SPGR(file,DicomFolder)


if ~isempty(strfind(DicomFolder,'/'))
    DicomFolderName = strfind(DicomFolder,'/');
    DicomFolderName = DicomFolder(DicomFolderName(end-1)+1:DicomFolderName(end)-1);
    Dicom_file_temp = dir([DicomFolder,'*.dcm']);
    Dicom_file_temp = [Dicom_file_temp(1).folder,'/',Dicom_file_temp(1).name];
else
    DicomFolderName = DicomFolder;
    Dicom_file_temp = DicomFolder;
end

load([file.folder,'/',file.name],'im_bin','para')
im_bin = crop_half_FOV(abs(im_bin));

orintation = orintation_detection(sum(im_bin(:,:,:,1),3));

im_bin = orintate_image(im_bin,orintation);


[sx,sy,nof,nm,ns] = size(im_bin);
im_bin = im_bin(:,:,:,:);
order = 1:nm:nm*ns;
for i=ns:-1:2
    order = [order;i:nm:nm*ns];
end
im_bin = im_bin(:,:,:,order);
Nphase = size(para.Recon.bins,1);
Ncycle = nof/Nphase;
im_bin = reshape(im_bin,[sx,sy,Nphase,Ncycle,nm*ns]);


%im_bin = permute(im_bin,[1,2,3,5,4]);
%im_bin = reshape(im_bin,sx,sy,nof,nm*ns);
im_bin = im_bin/max(im_bin(:))*600;%scale

seriesNum = file.name(4:8);
seriesDesc = para.dir.load_kSpace_name(1:end-14);
%seriesDesc = [seriesDesc,'F1000_T',num2str(para.weight_tTV*1000),'_S',num2str(para.weight_sTV*1000)];

parfor nSlice=1:nm*ns
    DicomFolderLocalSys = ['ReconData/',DicomFolderName,'_sys_slice_',num2str(nSlice,'%02.f')];
    mkdir(DicomFolderLocalSys)
    DicomFolderLocalDia = ['ReconData/',DicomFolderName,'_dia_slice_',num2str(nSlice,'%02.f')];
    mkdir(DicomFolderLocalDia)
    Image_one_slice = im_bin(:,:,:,:,nSlice);
    Image_sys = squeeze(Image_one_slice(:,:,1,:));
    Image_dia = squeeze(Image_one_slice(:,:,round((Nphase+1)/2),:));
    seriesNum_sys = [seriesNum,num2str(nSlice,'%02.f')];
    seriesNum_dia = [seriesNum,num2str(nSlice+50,'%02.f')];
    for nFrame = 1:nof/Nphase
        outnamesys = [DicomFolderLocalSys,'/Img_',num2str(nFrame,'%03.f'),'.dcm'];
        outnamedia = [DicomFolderLocalDia,'/Img_',num2str(nFrame,'%03.f'),'.dcm'];
        addDicomHeader(Dicom_file_temp, Image_sys(:,:,nFrame), seriesNum_sys, [seriesDesc,'sys'], nFrame, outnamesys)
        addDicomHeader(Dicom_file_temp, Image_dia(:,:,nFrame), seriesNum_dia, [seriesDesc,'dia'], nFrame, outnamedia)
    end
end

end