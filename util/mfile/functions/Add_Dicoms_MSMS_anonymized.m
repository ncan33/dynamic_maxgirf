function Add_Dicoms_MSMS_anonymized(file,DicomFolder)


DicomFolderName = DicomFolder;
Dicom_file_temp = DicomFolder;

load([file.folder,'/',file.name])

nof = size(Image,3);

if isfield(para{1}.Recon,'PD_frames')
    PD = para{1}.Recon.PD_frames;
    PD = find(PD);
    PD = PD(end);
else
    PD = 10;
end

Image(:,:,1:PD,:,:) = [];
if nof>80
    Image(:,:,end-19:end,:,:) = [];
end
[sx,sy,nof,nm,ns] = size(Image);

Image = permute(Image,[1,2,3,5,4]);
Image = reshape(Image,sx,sy,nof,nm*ns);
Image = Image/max(Image(:))*2000;%scale

seriesNum = file.name(4:8);
seriesDesc = [];

parfor nSlice=1:nm*ns
    DicomFolderLocal = ['/v/raid1b/ytian/MRIdata/TestData/122017_short_long_axis_Paper/Dicoms_All_new/',DicomFolderName,'_slice_',num2str(nSlice,'%02.f')];
    mkdir(DicomFolderLocal)
    Image_one_slice = Image(:,:,:,nSlice);
    seriesNum_temp = [seriesNum,num2str(nSlice,'%02.f')];
    for nFrame = 1:nof
        outname = [DicomFolderLocal,'/Img_',num2str(nFrame,'%03.f'),'.dcm'];
        data_in = Image_one_slice(:,:,nFrame);
        addDicomHeader(Dicom_file_temp, data_in, seriesNum_temp, seriesDesc, nFrame, outname)
    end
end

end