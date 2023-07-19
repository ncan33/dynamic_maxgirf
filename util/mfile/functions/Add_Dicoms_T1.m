function Add_Dicoms_T1(file,DicomFolder,Folder)

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
% T1 = rot90(T1,2);
[sx,sy,ns] = size(T1);

seriesNum = ['1',file.name(9:13)];
seriesDesc = ['MID',file.name(9:13)];
%seriesDesc = [seriesDesc,'F1000_T',num2str(para{1}.weight_tTV*1000),'_S',num2str(para{1}.weight_sTV*1000)];

for nSlice=1:ns
    DicomFolderLocal = ['ReconData/',Folder,'_',DicomFolderName,'_T1_slice_',num2str(nSlice,'%02.f')];
    mkdir(DicomFolderLocal)
    T1_one_slice = T1(:,:,nSlice);
    seriesNum_temp = [seriesNum,num2str(nSlice,'%02.f')];
    
    outname = [DicomFolderLocal,'/Img_',num2str(1,'%03.f'),'.dcm'];

    addDicomHeader(Dicom_file_temp, T1_one_slice, seriesNum_temp, seriesDesc, nSlice, outname)

end
