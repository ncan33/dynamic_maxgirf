function Add_Dicoms_3D(file,DicomFolder,name_add,seriesDesc,orintation)
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
        Dicom_file_temp = dir([DicomFolder,'IM*.dcm']);
        Dicom_file_temp = [Dicom_file_temp(1).folder,'/',Dicom_file_temp(1).name];
end

load([file.folder,'/',file.name])
if isfield(Image,'inj1')
    Image = cat(4,Image.pd,Image.inj1);
elseif isfield(Image,'inj2')
    Image = cat(4,Image.pd,Image.inj1,Image.inj2);
end
[sx,sy,ns,nof] = size(Image);
if ns>nof
    Image = permute(Image,[1,2,4,3]);
    [sx,sy,ns,nof] = size(Image);
end

Image = permute(Image,[1,2,4,3]);
Image = Image/max(Image(:))*4000;
% Image = orintate_image(Image,orintation);

seriesNum = file.name(4:8);
%seriesDesc = para.dir.load_kSpace_name(1:end-10);

seriesDesc = [seriesDesc,'F1000_T',num2str(para.weight_tTV*1000),'_S',num2str(para.weight_sTV*1000)];

parfor nSlice=1:ns
    DicomFolderLocal = ['ReconData/',DicomFolderName,name_add,'slice_',num2str(nSlice,'%02.f')];
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