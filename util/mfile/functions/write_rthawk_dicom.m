function write_rthawk_dicom(recon_file, dicom_dir)

%% recon dir
% recon_file = dir('./usc_disc_yt_2022_07_08_142402_multi_slice_golden_angle_spiral_ssfp_slice_30_fov_240_n63_tread3.mat');
write_name = recon_file.name(1:end-4);
recon_file = fullfile(recon_file.folder, recon_file.name);


%% dicom dir
% dicom_dir = './dicom_hawk/';

%% write dir
write_dir = sprintf('./dicom_recon/%s/', write_name);
if ~isfolder(write_dir)
    mkdir(write_dir);
end

%% find time of data
idx = strfind(recon_file, 'usc_disc_yt_'); 
time = recon_file(idx + 12 : idx + 28);

%% find how many subfolders in dicom
path_search_all = {};
path_search_all{1} = dicom_dir;
sub_folder_temp = search_sub_folder(path_search_all{1});
sub_folder_to_be_searched = sub_folder_temp;
sub_folder_all = [dicom_dir; sub_folder_temp];
while ~isempty(sub_folder_to_be_searched)
    sub_folder_temp = search_sub_folder(sub_folder_to_be_searched{1});
    sub_folder_to_be_searched(1) = [];
    sub_folder_to_be_searched = [sub_folder_to_be_searched; sub_folder_temp];
    sub_folder_all = [sub_folder_all; sub_folder_temp];
end

%% load first dicom from each sub folder
all_dicom = [];
for K = 1 : length(sub_folder_all)
    thisdir = sub_folder_all{K};
    dicom_temp = dir(fullfile(thisdir, '*.dcm'));
    if ~isempty(dicom_temp)
        all_dicom = [all_dicom; dicom_temp(1)];
    end
end

%% find time of each dicom
n_dicom = length(all_dicom);
for i = 1:n_dicom
    dicom_file = fullfile(all_dicom(i).folder, all_dicom(i).name);
    dinfo = dicominfo(dicom_file);
    dicom_time(i) = str2double(dinfo.SeriesTime(1:6));
end

%% find the best match
d = str2double(time(12:end)) - dicom_time;
[~, idx] = min(abs(d));

%% load the templete dicom info
dicom_file = fullfile(all_dicom(idx).folder, all_dicom(idx).name);
dinfo = dicominfo(dicom_file);

%% correct some time
load(recon_file)
nslice = size(image_recon, 3);
SeriesTime = str2double(dinfo.SeriesTime);

%% correct some slice location
q0 = para.kspace_info.user_QuaternionW;
q1 = para.kspace_info.user_QuaternionX;
q2 = para.kspace_info.user_QuaternionY;
q3 = para.kspace_info.user_QuaternionZ;
rot = [2 * (q0^2 + q1^2 ) - 1,   2 * (q1*q2 - q0*q3),    2 * (q1*q3 + q0*q2);
       2 * (q1*q2 + q0*q3),     2 * (q0^2 + q2^2 ) - 1,  2 * (q2*q3 - q0*q1);
       2 * (q1*q3 - q0*q2),     2 * (q2*q3 + q0*q1),    2 * (q0^2 + q3^2) - 1];
   
x0 = dinfo.ImagePositionPatient(1);
y0 = dinfo.ImagePositionPatient(2);
z0 = dinfo.ImagePositionPatient(3);

SliceGap = para.kspace_info.user_SliceGap;

SliceShifts = ((1:nslice) - floor(nslice / 2) - 1) * SliceGap;
% SliceShifts = ((1:nslice) - 1) * SliceGap;
SliceShiftsRot = rot * [zeros(2, nslice); SliceShifts];

%% write dicom
image_recon = int16(image_recon * 32767 / 40);
time_per_slice = para.kspace_info.user_TR * para.Recon.narm / 1000 / 1000; % [sec]

SliceLocation0          = dinfo.SliceLocation;
ImagePositionPatient0   = dinfo.ImagePositionPatient;

for i = 1:nslice
    warning off
    AcquisitionTime         = num2str(round((SeriesTime + time_per_slice * i) * 1000) / 1000);
    ContentTime             = AcquisitionTime;
    ImagesInAcquisition     = nslice;
    InstanceNumber          = i;
    FrameTime               = time_per_slice * 1000;
    SliceLocation           = SliceLocation0 + SliceShifts(i);
    ImagePositionPatient    = ImagePositionPatient0 + SliceShiftsRot(:, i);
    
    
    dinfo.AcquisitionTime       = AcquisitionTime;
    dinfo.ContentTime           = ContentTime;
    dinfo.ImagesInAcquisition   = ImagesInAcquisition;
    dinfo.InstanceNumber        = InstanceNumber;
    dinfo.FrameTime             = FrameTime;
    dinfo.SliceLocation         = SliceLocation; 
    dinfo.ImagePositionPatient  = ImagePositionPatient;
    
    file_name = sprintf('%s%s_slice_%03g.dcm', write_dir, write_name, i);
    
    dicomwrite(image_recon(:, :, i), file_name, dinfo);
    
end

end

