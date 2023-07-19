function Add_Dicoms_MSMS_self_gating(file,DicomFolder,ShortLongFlag,Folder)

if nargin == 2
    ShortLongFlag = 0;
    Folder = [];
elseif nargin == 3
    Folder = [];
end

DicomFolderName = strfind(DicomFolder,'/');
DicomFolderName = DicomFolder(DicomFolderName(end-1)+1:DicomFolderName(end)-1);

load([file.folder,'/',file.name])
[sx,sy,nof,nset,nsms] = size(Image);

if ShortLongFlag
    Image(:,:,:,ShortLongFlag,:) = fliplr(Image(:,:,:,ShortLongFlag,:));
end
%Image = permute(Image,[1,2,3,5,4]);
%Image = reshape(Image,sx,sy,nof,nm*ns);
Image = Image/max(Image(:))*2000;%scale

Dicom_file_temp = dir([DicomFolder,'*.dcm']);
Dicom_file_temp = [Dicom_file_temp(1).folder,'/',Dicom_file_temp(1).name];

seriesNum = file.name(4:8);
seriesDesc = para{1}.dir.load_kSpace_name(1:end-10);
seriesDesc = [seriesDesc,'F1000_T',num2str(para{1}.weight_tTV*1000),'_S',num2str(para{1}.weight_sTV*1000)];

for i=1:nset
    sys(i,:) = logical(sum(para{i}.Recon.bins(1:3:end,:),1));
    dia(i,:) = logical(sum(para{i}.Recon.bins(2:3:end,:),1));
end
Nf = min(min(sum(sys,2)),min(sum(dia,2)));
for i=1:nset
    temp = find(sys(i,:));
    sys_short(i,:) = temp(1:Nf);
    temp = find(dia(i,:));
    dia_short(i,:) = temp(1:Nf);
end
Image(:,:,1:10,:,:) = [];
for i=1:nset
    Image_short(:,:,:,i,:,1) = Image(:,:,sys_short(i,:),i,:);
    Image_short(:,:,:,i,:,2) = Image(:,:,dia_short(i,:),i,:);
end
Image_short = permute(Image_short,[1,2,3,5,4,6]);
Image_short = reshape(Image_short,sx,sy,Nf,nset*nsms,2);


parfor nSlice=1:nset*nsms
    DicomFolderLocal_sys = ['ReconData/',Folder,'Self_Gating/sys/',DicomFolderName,'_slice_',num2str(nSlice,'%02.f')];
    mkdir(DicomFolderLocal_sys)
    Image_one_slice_sys = Image_short(:,:,:,nSlice,1);
    %seriesNum_temp = [seriesNum,num2str(nSlice,'%02.f')];
    
    DicomFolderLocal_dia = ['ReconData/',Folder,'Self_Gating/dia/',DicomFolderName,'_slice_',num2str(nSlice,'%02.f')];
    mkdir(DicomFolderLocal_dia)
    Image_one_slice_dia = Image_short(:,:,:,nSlice,2);
    seriesNum_temp = [seriesNum,num2str(nSlice,'%02.f')];
    
    for nFrame = 1:Nf
        outname_sys = [DicomFolderLocal_sys,'/Img_',num2str(nFrame,'%03.f'),'.dcm'];
        data_in_sys = Image_one_slice_sys(:,:,nFrame);
        addDicomHeader(Dicom_file_temp, data_in_sys, seriesNum_temp, seriesDesc, nFrame, outname_sys)
        
        outname_dia = [DicomFolderLocal_dia,'/Img_',num2str(nFrame,'%03.f'),'.dcm'];
        data_in_dia = Image_one_slice_dia(:,:,nFrame);
        addDicomHeader(Dicom_file_temp, data_in_dia, seriesNum_temp, seriesDesc, nFrame, outname_dia)
    end
end

end