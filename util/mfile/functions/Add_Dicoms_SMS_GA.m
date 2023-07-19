function Add_Dicoms_SMS_GA(file,DicomFolder,name_add,seriesDesc)
switch DicomFolder
    case 'no'
        Dicom_file_temp = 'no';
        DicomFolderName = file.name;
        idx = strfind(DicomFolderName,'MID');
        DicomFolderName = DicomFolderName(idx:end-18);
        idx = strfind(file.folder,'/');
        Dicom_file_temp = [file.folder(idx(7)+1:idx(8)-1),'_',DicomFolderName];
    otherwise
        DicomFolderName = file.name;
        idx = strfind(DicomFolderName,'MID');
        DicomFolderName = DicomFolderName(idx:end-18);
        %DicomFolderName = strfind(DicomFolder,'/');
        %DicomFolderName = DicomFolder(DicomFolderName(end-1)+1:DicomFolderName(end)-1);
        Dicom_file_temp = dir([DicomFolder,'*.dcm']);
        Dicom_file_temp = [Dicom_file_temp(1).folder,'/',Dicom_file_temp(1).name];
end

load([file.folder,'/',file.name])
[sx,sy,nof,nm,ns] = size(Image);

Image = permute(Image,[1,2,3,5,4]);
Image = reshape(Image,sx,sy,nof,nm*ns);

seriesNum = file.name(4:8);
% if length(para) > 1
%     para = para{1};
% end
% %seriesDesc = para.dir.load_kSpace_name(1:end-10);
% seriesDesc = [seriesDesc,'F1000_T',num2str(para.weight_tTV*1000),'_S',num2str(para.weight_sTV*1000)];

parfor nSlice=1:nm*ns
    DicomFolderLocal = ['ReconData/',DicomFolderName,name_add,'slice_',num2str(nSlice,'%02.f')];
    mkdir(DicomFolderLocal)
    Image_one_slice = Image(:,:,:,nSlice);
    Image_one_slice = Image_one_slice/max(Image_one_slice(:))*4000;
    seriesNum_temp = [seriesNum,num2str(nSlice,'%02.f')];
    for nFrame = 1:nof
        outname = [DicomFolderLocal,'/Img_',num2str(nFrame,'%03.f'),'.dcm'];
        data_in = Image_one_slice(:,:,nFrame);
        addDicomHeader(Dicom_file_temp, data_in, seriesNum_temp, seriesDesc, nFrame, outname)
    end
end
end