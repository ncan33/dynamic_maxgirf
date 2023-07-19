%ADDDICOMHEADER
%   Function to add meta data from a given file, to a new dicom file that
%   is being created, with the specified tags modified.
%
%   Format: ADDDICOMHEADER(templateFilename, data_in,seriesNum, outname)
%       where, templateFilename is the name of the file whose dicom tags
%       have to be copied
%       data_in is the input data that goes into the new dicom file thats
%       being written
%       seriesNum is the value of the series number tag that needs to be
%       replaced from the original dicom tag
%       seriesDesc is the series description that needs to be replaced from
%       the original dicom tag
%       instanceNum is the value of instance number that needs to be used
%       in the new dicom file
%       outname is the output file name of the new dicom file being
%       created.
function addDicomHeader(templateFilename, data_in, seriesNum, seriesDesc, instanceNum, outname)

%   Copyright UCAIR, University of Utah, 2007.
%   Author: Edward DiBella
%   Last Modified by: Sathya Vijayakumar Dec 6th 2007.

if isempty(dir(templateFilename))
    md = create_fake_dicom_info(templateFilename);
else
    md = dicominfo(templateFilename);
end
if exist('seriesDesc')
    md.SeriesDescription = seriesDesc;
end
md.SeriesNumber = str2num(seriesNum);

data = data_in;
flatdata = reshape(data,1,numel(data));

md.Height=size(data,1);
md.Width=size(data,2);
outdata=double(data);

range=[];
if (isempty(range))
    range = [min(flatdata) max(flatdata)];
    scale = min(65535,max(flatdata));
end
badData = find(~isfinite(outdata));

%scale= 65535/max(flatdata);
% scale= 511/max(flatdata);  % more like Siemens range, trying this 7/15/09 EVRD to see if quantization noise
%scale= max_int/max(flatdata);
%scale=1;

%outdata = uint16(outdata/scale*700);
outdata = uint16(outdata);
outdata(badData) = range(1);

md.BitDepth = 16;
md.BitsAllocated = 16;
md.BitsStored = 16;

md.SmallestImagePixelValue = scale*min(flatdata);
md.LargestImagePixelValue = scale*max(flatdata);
md.InstanceNumber = instanceNum;

%if strfind(templateFilename,'stress')
%    md.StudyDescription = 'stress';
%elseif strfind(templateFilename,'rest')
%    md.StudyDescription = 'rest';
%end
warning off
dicomwrite(outdata,outname,md);    % was outdata'
