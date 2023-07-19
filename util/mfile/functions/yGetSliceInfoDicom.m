function SliceInfo = yGetSliceInfoDicom(DicomFolder,nSET)

if DicomFolder(end) ~= '/'
    DicomFolder(end+1) = '/';
end

a = dir([DicomFolder,'*.dcm']);

SliceInfo.PixelSpacing = zeros(length(a),2,'single');
SliceInfo.ImageOrientation = zeros(length(a),6,'single');
SliceInfo.ImagePosition = zeros(length(a),3,'single');
SliceInfo.SliceThickness = zeros(length(a),1,'single');
SliceInfo.SlicePosition = zeros(length(a),1,'single');
SliceInfo.SpacingBetweenSlices = zeros(length(a),1,'single');
SliceInfo.FOV = zeros(length(a),2,'single');

for i=1:nSET
    dinfo_temp = dicominfo([a(i).folder,'/',a(i).name]);
    SliceInfo.PixelSpacing(i,:) = dinfo_temp.PixelSpacing;
    SliceInfo.ImageOrientation(i,:) = dinfo_temp.ImageOrientationPatient;
    SliceInfo.ImagePosition(i,:) = dinfo_temp.ImagePositionPatient;
    SliceInfo.SliceThickness(i,:) = dinfo_temp.SliceThickness;
    SliceInfo.SlicePosition(i,:) = dinfo_temp.SliceLocation;
    SliceInfo.SpacingBetweenSlices(i,:) = dinfo_temp.SpacingBetweenSlices;
    SliceInfo.FOV(i,1) = str2num(dinfo_temp.Private_0051_100c(5:7));
    SliceInfo.FOV(i,2) = str2num(dinfo_temp.Private_0051_100c(9:11));
end