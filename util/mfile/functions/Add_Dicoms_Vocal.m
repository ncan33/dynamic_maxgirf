function Add_Dicoms_Vocal(file,DicomFolder,ShortLongFlag)

if nargin==2
    ShortLongFlag = false;
end

DicomFolderName = strfind(DicomFolder,'/');
DicomFolderName = DicomFolder(DicomFolderName(end-1)+1:DicomFolderName(end)-1);

load([file.folder,'/',file.name])
[sx,sy,nof,nm,ns] = size(Image);
if ShortLongFlag
    Image(:,:,:,1,:) = rot90(Image(:,:,:,1,:),2);
end
Image = reshape(Image,sx,sy,nof,nm*ns);

%MID = strfind(file.name,'MID');
%MID = file.name(MID+3:MID+5);

%raid1apath = pwd;
%s = strfind(raid1apath,'ytian');
%raid1apath = [raid1apath(1:s-1),'gadluru',raid1apath(s+5:end),'/'];

Dicom_file_temp = dir([DicomFolder,'*.dcm']);
Dicom_file_temp = [Dicom_file_temp(1).folder,'/',Dicom_file_temp(1).name];

seriesNum = file.name(4:8);
seriesDesc = para{1}.dir.load_kSpace_name(1:end-10);
seriesDesc = [seriesDesc,'F1000_T',num2str(para{1}.weight_tTV*1000),'_S',num2str(para{1}.weight_sTV*1000)];

%outname_temp = [pwd,'ReconData'];

%DicomFolderLocal = [DicomFolderName,'_slice_',num2str(nSlice,'%02.f')];
DicomFolderLocal = ['ReconData/',DicomFolderName];
mkdir(DicomFolderLocal)
seriesNum_temp = seriesNum;
for nFrame = 1:nof
    outname = [DicomFolderLocal,'/Img_',num2str(nFrame,'%03.f'),'.dcm'];
    data_in = Image(:,:,nFrame);
    addDicomHeader(Dicom_file_temp, data_in, seriesNum_temp, seriesDesc, nFrame, outname)
end

end