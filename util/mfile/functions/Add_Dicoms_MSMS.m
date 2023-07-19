function Add_Dicoms_MSMS(file,DicomFolder,ShortLongFlag,Folder)

if nargin == 2
    ShortLongFlag = 0;
    Folder = [];
elseif nargin == 3
    Folder = [];
end

if ~isempty(strfind(DicomFolder,'/'))
    DicomFolderName = strfind(DicomFolder,'/');
    DicomFolderName = DicomFolder(DicomFolderName(end-1)+1:DicomFolderName(end)-1);
    Dicom_file_temp = dir([DicomFolder,'*.dcm']);
    Dicom_file_temp = [Dicom_file_temp(1).folder,'/',Dicom_file_temp(1).name];
else
    DicomFolderName = DicomFolder;
    Dicom_file_temp = DicomFolder;
end

load([file.folder,'/',file.name])
Image = rot90(Image);
[sx,sy,nof,nm,ns] = size(Image);
if ShortLongFlag
    Image(:,:,:,ShortLongFlag,:) = fliplr(Image(:,:,:,ShortLongFlag,:));
end
Image = permute(Image,[1,2,3,5,4]);
Image = reshape(Image,sx,sy,nof,nm*ns);
Image = Image/max(Image(:))*4000;%scale

seriesNum = file.name(4:8);
seriesDesc = para{1}.dir.load_kSpace_name(1:end-10);
seriesDesc = [seriesDesc,'F1000_T',num2str(para{1}.weight_tTV*1000),'_S',num2str(para{1}.weight_sTV*1000)];

parfor nSlice=1:nm*ns
    DicomFolderLocal = ['ReconData/',Folder,DicomFolderName,'_slice_',num2str(nSlice,'%02.f')];
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