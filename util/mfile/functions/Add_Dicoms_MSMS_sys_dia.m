function Add_Dicoms_MSMS_sys_dia(file,DicomFolder,ShortLongFlag,Folder)

if nargin == 2
    ShortLongFlag = 0;
    Folder = [];
elseif nargin == 3
    Folder = [];
end

switch DicomFolder
    case 'no'
        Dicom_file_temp = 'no';
        DicomFolderName = file.name;
        idx = strfind(DicomFolderName,'MID');
        DicomFolderName = DicomFolderName(idx:end-18);
        idx = strfind(file.folder,'/');
        Dicom_file_temp = [file.folder(idx(7)+1:idx(8)-1),'_',DicomFolderName];
    otherwise
        %DicomFolderName = strfind(DicomFolder,'/');
        %DicomFolderName = DicomFolder(DicomFolderName(end-1)+1:DicomFolderName(end)-1);
        DicomFolderName = file.name;
        idx = strfind(DicomFolderName,'MID');
        DicomFolderName = DicomFolderName(idx:end-18);
        Dicom_file_temp = dir([DicomFolder,'*.dcm']);
        Dicom_file_temp = [Dicom_file_temp(1).folder,'/',Dicom_file_temp(1).name];
end
%sys
load([file.folder,'/',file.name])
if length(para) == 1
    para_temp = para;
    clear para
    para{1} = para_temp;
end
[sx,sy,nof,nm,ns] = size(Image_sys);

orintation = orintation_detection(sum(Image_sys(:,:,:,1),3));
Image_sys = orintate_image(Image_sys,orintation);
Image_dia = orintate_image(Image_dia,orintation);

if ShortLongFlag
    Image_sys(:,:,:,ShortLongFlag,:) = fliplr(Image_sys(:,:,:,ShortLongFlag,:));
end
Image_sys = permute(Image_sys,[1,2,3,5,4]);
Image_sys = reshape(Image_sys,sx,sy,nof,nm*ns);
Image_sys = Image_sys/max(Image_sys(:))*200;%scale

seriesNum = ['1',file.name(4:8)];
seriesDesc = ['MID',para{1}.dir.load_kSpace_name(11:13),para{1}.dir.load_kSpace_name(23:end-10)];
seriesDesc = [seriesDesc,'F1000_T',num2str(para{1}.weight_tTV*1000),'_S',num2str(para{1}.weight_sTV*1000)];

parfor nSlice=1:nm*ns
    DicomFolderLocal = ['ReconData/',Folder,DicomFolderName,'_sys_slice_',num2str(nSlice,'%02.f')];
    mkdir(DicomFolderLocal)
    Image_sys_one_slice = Image_sys(:,:,:,nSlice);
    seriesNum_temp = [seriesNum,num2str(nSlice,'%02.f')];
    for nFrame = 1:nof
        outname = [DicomFolderLocal,'/Img_',num2str(nFrame,'%03.f'),'.dcm'];
        data_in = Image_sys_one_slice(:,:,nFrame);
        addDicomHeader(Dicom_file_temp, data_in, seriesNum_temp, seriesDesc, nFrame, outname)
    end
end
%dia
[sx,sy,nof,nm,ns] = size(Image_dia);
if ShortLongFlag
    Image_dia(:,:,:,ShortLongFlag,:) = fliplr(Image_dia(:,:,:,ShortLongFlag,:));
end
Image_dia = permute(Image_dia,[1,2,3,5,4]);
Image_dia = reshape(Image_dia,sx,sy,nof,nm*ns);
Image_dia = Image_dia/max(Image_dia(:))*200;%scale

seriesNum = ['2',file.name(4:8)];
seriesDesc = ['MID',para{1}.dir.load_kSpace_name(11:13),para{1}.dir.load_kSpace_name(23:end-10)];
seriesDesc = [seriesDesc,'F1000_T',num2str(para{1}.weight_tTV*1000),'_S',num2str(para{1}.weight_sTV*1000)];

parfor nSlice=1:nm*ns
    DicomFolderLocal = ['ReconData/',Folder,DicomFolderName,'_dia_slice_',num2str(nSlice,'%02.f')];
    mkdir(DicomFolderLocal)
    Image_dia_one_slice = Image_dia(:,:,:,nSlice);
    seriesNum_temp = [seriesNum,num2str(nSlice,'%02.f')];
    for nFrame = 1:nof
        outname = [DicomFolderLocal,'/Img_',num2str(nFrame,'%03.f'),'.dcm'];
        data_in = Image_dia_one_slice(:,:,nFrame);
        addDicomHeader(Dicom_file_temp, data_in, seriesNum_temp, seriesDesc, nFrame, outname)
    end
end

end