function Image_reg = ANTS_2D(Image_moving,Image_ref)
%
% ccc
%
% load('ReconData/SET13692_meas_MID00075_FID26371_UCAIR_Ungated2D_stress_11ml.mat','Image_sys')
%
% test1 = double(Image_sys(:,:,1,1));
% test2 = double(Image_sys(:,:,1,2));
[sx,sy,nof] = size(Image_moving);
Image_reg = zeros([sx,sy,nof]);
mkdir ANTS_temp/
working_dir = [pwd,'/ANTS_temp/'];

for i=1:nof
%     [sx,sy] = size(test1);
    
    filename_ref_raw = strcat(working_dir,'/ref_',int2str(i),'.raw');
    fid = fopen (filename_ref_raw, 'wb', 'ieee-le');
    fwrite (fid, Image_moving(:,:,i), 'float');
    fclose(fid);
    
    filename_ref = strcat(working_dir,'/ref_',int2str(i),'.mhd');
    fid = fopen (filename_ref, 'wb');
    header_info = sprintf('ObjectType = Image\nNDims = 2\nBinaryData = True\nBinaryDataByteOrderMSB = False\nDimSize = %d %d\nElementType = MET_FLOAT\nElementDataFile = %s\n',sx, sy, filename_ref_raw);
    fprintf(fid,'%s',header_info);
    fclose(fid);
    
    filename_tar_raw = strcat(working_dir,'/tar_',int2str(i),'.raw');
    fid = fopen (filename_tar_raw, 'wb', 'ieee-le');
    fwrite (fid, Image_ref(:,:,i), 'float');
    fclose(fid);
    
    filename_tar = strcat(working_dir,'/tar_',int2str(i),'.mhd');
    fid = fopen (filename_tar, 'wb');
    header_info = sprintf('ObjectType = Image\nNDims = 2\nBinaryData = True\nBinaryDataByteOrderMSB = False\nDimSize = %d %d\nElementType = MET_FLOAT\nElementDataFile = %s\n',sx, sy, filename_tar_raw);
    fprintf(fid,'%s',header_info);
    fclose(fid);
    
    filename_out = strcat(working_dir,'/out_',int2str(i));
    
    % don't know what are these yet
    step = 0.1;
    a1 = 3;
    a2 = 0;
    
    eval_command = horzcat('!/v/raid1b/gadluru/softs/ANTS/092514/antsbin/bin/antsRegistration -d 2 -r \[',filename_ref,' , ',filename_tar,', 1\]',...
        ' -m mattes\[',filename_ref,',',filename_tar,', 1, 32, regular, 0.1\]',...
        ' -t affine\[0.1\]',...
        ' -c \[100x100x10, 1e-30, 4\]',...
        ' -s 4x2x1vox',...
        ' -f 3x2x1',...
        ' -u -histogram-matching 1',...
        ' -m demons\[',filename_ref,' , ',filename_tar,', 1, 1\]',...
        [' -t SyN\[',num2str(step),',',num2str(a1) ,', ',num2str(a2),'\]'],...
        ' -c \[100x70x50, 1e-8, 5\]',...
        ' -s 1x0.5x0vox',...
        ' -f 1x1x1',...
        ' -o ',filename_out);
%     eval_command = horzcat('!/v/raid1b/gadluru/softs/ANTS/092514/antsbin/bin/antsRegistration -d 2 -r \[',filename_ref,' , ',filename_tar,', 1\]',...
%         ' -u -histogram-matching 1',...
%         ' -m CC\[',filename_ref,' , ',filename_tar,', 1, 4\]',...
%         [' -t SyN\[',num2str(step),',',num2str(a1) ,', ',num2str(a2),'\]'],...
%         ' -c \[100x70x50, 1e-8, 5\]',...
%         ' -s 1x0.5x0vox',...
%         ' -f 2x1x1',...
%         ' -o ',filename_out);
    eval(eval_command)
    
    filename_out_warp=strcat(filename_out,'1Warp.nii.gz');
    filename_out_affine=strcat(filename_out,'0GenericAffine.mat');
    filename_tar_warp_mhd=strcat(filename_tar(1:end-4),'_warp.mhd');
    
    eval_command2 = horzcat('!/v/raid1b/gadluru/softs/ANTS/092514/antsbin/bin/antsApplyTransforms -d 2 -i ',filename_tar,' -o ',filename_tar_warp_mhd,' -r ',filename_ref,'  -n BSpline -t ',filename_out_affine,' -t ',filename_out_warp);
    eval(eval_command2)

    filename_tar_warp_raw=strcat(filename_tar(1:end-4),'_warp.raw');
    fid = fopen (filename_tar_warp_raw, 'r', 'ieee-le');
    
    temp_var=fread (fid,'double');
    
    fclose (fid);
    
    
    Image_reg(:,:,i)=reshape(temp_var,[sx sy]);
end