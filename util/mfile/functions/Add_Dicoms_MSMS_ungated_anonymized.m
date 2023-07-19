function Add_Dicoms_MSMS_ungated_anonymized(file,DicomFolder)


DicomFolderName = DicomFolder;
Dicom_file_temp = DicomFolder;

load([file.folder,'/',file.name])

% if isfield(para{1}.Recon,'PD_frames')
%     PD = para{1}.Recon.PD_frames;
%     PD = find(PD);
%     PD = PD(end);
% else
%     PD = 10;
% end

% cardiac_signal(end-19:end,:) = [];
% Image(:,:,1:PD,:,:) = [];

[sx,sy,nof_sys,nset,nsms] = size(Image_sys);
Image_sys = Image_sys/max(Image_sys(:))*2000;%scale


[sx,sy,nof_dia,nset,nsms] = size(Image_dia);
Image_sys = Image_sys/max(Image_sys(:))*2000;%scale

nof = min(nof_sys,nof_dia);

Image_sys(:,:,nof-19:end,:,:) = [];
Image_dia(:,:,nof-19:end,:,:) = [];
nof = nof-20;

Image_sys = permute(Image_sys,[1,2,3,5,4]);
Image_sys = reshape(Image_sys,[sx,sy,nof,nset*nsms]);

Image_dia = permute(Image_dia,[1,2,3,5,4]);
Image_dia = reshape(Image_dia,[sx,sy,nof,nset*nsms]);
%Image = permute(Image,[1,2,3,5,4]);
%Image = reshape(Image,sx,sy,nof,nsms*nset);

seriesNum = file.name(4:8);
seriesDesc = [];
% fprintf('self-gating: \n')
% for i=1:nset
%     sys(i,:) = cardiac_signal(:,i)==1;
%     dia(i,:) = cardiac_signal(:,i)==2;
%     fprintf(['set ',num2str(i,'%02.f'),' ',num2str(100-(sum(sys(i,:)) + sum(dia(i,:)))./nof*100,'%02.1f'),'%% dropped \n']);
% end

% Nf = min(min(sum(sys,2)),min(sum(dia,2)));
% for i=1:nset
%     temp = find(sys(i,:));
%     sys_short(i,:) = temp(1:Nf);
%     temp = find(dia(i,:));
%     dia_short(i,:) = temp(1:Nf);
% end
% 
% for i=1:nset
%     Image_short(:,:,:,i,:,1) = Image(:,:,sys_short(i,:),i,:);
%     Image_short(:,:,:,i,:,2) = Image(:,:,dia_short(i,:),i,:);
% end
% Image_short = permute(Image_short,[1,2,3,5,4,6]);
% Image_short = reshape(Image_short,sx,sy,Nf,nset*nsms,2);

parfor nSlice=1:nset*nsms
    DicomFolderLocal_sys = ['/v/raid1b/ytian/MRIdata/TestData/122017_short_long_axis_Paper/Dicoms_All_new/',DicomFolderName,'_sys_slice_',num2str(nSlice,'%02.f')];
    mkdir(DicomFolderLocal_sys)
    Image_one_slice_sys = Image_sys(:,:,:,nSlice);
    seriesNum_temp_sys = [seriesNum,num2str(nSlice,'%02.f')];
    
    DicomFolderLocal_dia = ['/v/raid1b/ytian/MRIdata/TestData/122017_short_long_axis_Paper/Dicoms_All_new/',DicomFolderName,'_dia_slice_',num2str(nSlice,'%02.f')];
    mkdir(DicomFolderLocal_dia)
    Image_one_slice_dia = Image_dia(:,:,:,nSlice);
    seriesNum_temp_dia = [seriesNum,num2str(nSlice+nset*nsms,'%02.f')];
    
    for nFrame = 1:nof
        outname_sys = [DicomFolderLocal_sys,'/Img_',num2str(nFrame,'%03.f'),'.dcm'];
        data_in_sys = Image_one_slice_sys(:,:,nFrame);
        addDicomHeader([Dicom_file_temp,' systole'], data_in_sys, seriesNum_temp_sys, seriesDesc, nFrame, outname_sys)
        
        outname_dia = [DicomFolderLocal_dia,'/Img_',num2str(nFrame,'%03.f'),'.dcm'];
        data_in_dia = Image_one_slice_dia(:,:,nFrame);
        addDicomHeader([Dicom_file_temp,' diastole'], data_in_dia, seriesNum_temp_dia, seriesDesc, nFrame, outname_dia)
    end
end

end