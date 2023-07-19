
% inputs - 
% template_dicom_directory='/v/raid1a/gadluru/MRIdata/Cardiac/Prisma/D090916_PD2S1/DicomData/CV_SMS_ungated_coverage_4ml_0_05dose_Real(122)';%strcat(template_dicom_directory,partname,int2str(tmpcounter),')');
% flag_orient

function Add_Dicom_header_3D(flag_orient,cardiac_phase,template_dicom_directory_path)

addpath /v/raid1a/gadluru/MRIdata/mat_files/

% findMIDrelations_CAIPI
% findMIDrelations_VD

% manual inputs - numcoils and MIDs

numCoil = 1:8;

filenames=dir('Pre_Interp*mat');
MIDs=zeros(length(filenames),1);
for i=1:length(filenames)
    a=filenames(i).name;
    b=strfind(a,'_');
    MIDs(i)=str2double(a(b(3)+1:b(4)-1));
end
mids=unique(MIDs);
mids(isnan(mids))=[];

% patient number
%tmp=pwd;
%tmp2=strfind(tmp,'ytian');
%raid1path= ['/v/raid1a/gadluru',tmp(tmp2+5:end),'/DicomData/'];

tmp = strfind(template_dicom_directory_path,'/');
part_seriesdesc = template_dicom_directory_path(tmp(end-2)+1:tmp(end-1)-1); % finding the current patient number

%{
try
    matlabpool close
catch err
    disp(err.message)
end

try
    matlabpool local 8
catch err
    disp(err.message)
end

try
    delete(gcp('nocreate'))
catch err
    disp(err.message)
end

try
    parpool('local',12)
catch err
    disp(err.message)
end
%}
%for gMID=1:length(mids) %662 663
    
    MID=mids;
    
    %     raid1path1=strcat(raid1path,'ReconData/mat_files/',part_seriesdesc,'_MID');
    
    %     filename2=strcat(raid1path1,int2str(MID),'_nufft_step',int2str(1000*stepsize),'_Coil1_slice1_iter',int2str(noi),'_F',int2str(1000*weight_fidelity),'_T',int2str(1000*weight_laplacian),'.mat')
    filename2=strcat('Pre_Interp_MID_',num2str(MID,'%03.f'),'_slice_',int2str(1),'_MC.mat');
    load(filename2);  % just to get size
    %img_est=recon_data_inter;
    
    %scale_factor=700/max(recon_data_inter(:)); % 06/08/17
    scale_factor=4000/max(recon_data_inter(:));
    %[sx sy sz]=size(img_est);
    
    %if (sz > 40)
    %    tframe=20;
    %else
    %    tframe=6;
    %end
    
    clear img_est
    
    %     tmp_filename=strcat(raid1path,'Processing/MIDinfolookup.mat');
    %     load(tmp_filename)
    
%     load MIDinfolookup.mat
    
    %current_mid=MID;%MIDinfolookup(MID);
    
    numSlice=1:length(dir(strcat('Pre_Interp_MID_',num2str(MID,'%03.f'),'_slice_*mat'))); %1:length(current_mid.SeriesNumbers); %%% this must be correctly specified for the data being combined

    for islice=numSlice
        try
        call_this(MID,islice,part_seriesdesc,scale_factor,flag_orient,cardiac_phase,template_dicom_directory_path)
        catch
        end
        
    end
%end % for MID
return;


function call_this(MID,islice,part_seriesdesc,scale_factor,flag_orient,cardiac_phase,template_dicom_directory_path)
weight_fidelity=1;
weight_laplacian=0.01; % regularization parameter 1
weight_TV=0.001; % regularization parameter 2
noi=150;

%sos=zeros(sx,sy,sz);
%nmcoil=0;%thk

fprintf(['slice number ',num2str(islice),'\n'])

filename=strcat('Pre_Interp_MID_',num2str(MID,'%03.f'),'_slice_',int2str(islice),'_MC.mat');
load(filename)

imgrr=recon_data_inter*scale_factor;
%imgrr=(sos/(max(max(max(sos)))))*500; %thk

%============================================thk_b


%             figure(70); colormap gray;
%             imagesc((abs(imgrr(:,:,tframe))));
%
%             fullimagetitle=strcat('Combined image','/ slice',int2str(islice),' / frame',int2str(tframe));
%             title(fullimagetitle);
%             axis image on;
%             colorbar('location','EastOutside');

%============================================thk_e

%%% write out dicoms instead of saving .mats ------ GA
tmpg=imgrr;%clear imgrr
imgrr=tmpg(:,:,2:end-3);

%template_dicom_directory=raid1path;


if strfind(cardiac_phase,'diastole')
    partname=['3D_perfusion_',cardiac_phase];%(current_mid.SeriesDescription);
elseif strfind(cardiac_phase,'systole')
    partname=['3D_perfusion_',cardiac_phase];%(current_mid.SeriesDescription);  
end

% partname_dia=strcat(current_mid.SeriesDescription,'_dia');
if(partname(end)~='(')
    partname=strcat(partname,'(');
%     partname_dia=strcat(partname_dia,'(');
end

tmpcounter=MID; % CAIPI % current_mid.SeriesNumbers(islice);

if strfind(cardiac_phase,'diastole')
tmpcounter2=tmpcounter*400+6;%+islice-1;
elseif strfind(cardiac_phase,'systole')
tmpcounter2=tmpcounter*300+6;%+islice-1;
end

%template_dicom_directory = template_dicom_directory_path;
%template_dicom_directory=strcat(template_dicom_directory,template_dicom_directory_path);

partname=strcat(partname,int2str(tmpcounter2+islice-1),')/');
% partname_dia=strcat(partname_dia,int2str(tmpcounter2+islice-1),')/');

seriesnum=tmpcounter2+islice-1;


%outname=strcat(raid1path,'ReconData/');
outname=strcat('./',part_seriesdesc,'/ReconData/');

outname=strcat(outname,partname);

mkdir(outname)

delete_previous=strcat(outname,'*dcm');
delete(delete_previous)

outname=strcat(outname,'MID_',num2str(MID,'%03.f'),'_');
outname=strcat(outname,'STCR_slice_',int2str(islice),'_iter',int2str(noi),'_F',int2str(1000*weight_fidelity),'_T',int2str(1000*weight_laplacian),'_S',int2str(1000*weight_TV),'_');
seriesdesc=strcat(part_seriesdesc,'_series',int2str(seriesnum),'_imgSTCR_slice',int2str(islice));

func_dicomHeader(template_dicom_directory_path, imgrr, outname, seriesnum, seriesdesc,flag_orient,cardiac_phase);
%%% write out dicoms instead of saving .mats done ------ GA

return;





function func_dicomHeader(templatedirtxt, indata, outname, seriesnum, seriesdesc,flag_orient,cardiac_phase);
%  func_dicomHeader(templatedirtxt, indata, outname, seriesnum, seriesdesc)
%   finds all files with suffix .dcm in given templatedir
% revised with more error checking  2/24/10

tmpp=strcat(templatedirtxt, '/*.dcm');
h = dir(tmpp);
slice_loc = [];
n=0;

for i = 1:1  %42 %-- Here the assumption is for our scans, one slice per folder -- %
    %h(i).name
    
    try
        d = dicominfo(strcat(templatedirtxt,'/',h(i).name));
        temp_num=d.ImageOrientationPatient;
        temp2_num=d.ImagePositionPatient;
        
        temp_a=length(find(temp_num<0));
        temp_b=length(find(temp2_num<0));
        
        %         if(temp_a==4&&temp_b==2 || temp_a==3&&temp_b==1)
        %             for temp_cntr=1:size(indata,3)
        %                 indata(:,:,temp_cntr)=(rot90(indata(:,:,temp_cntr),-1));
        %             end
        %         elseif(temp_a==2&&temp_b==2 || temp_a==3&&temp_b==2 || temp_a==4&&temp_b==1 || temp_a==2&&temp_b==3)
        %             for temp_cntr=1:size(indata,3)
        %                 indata(:,:,temp_cntr)=fliplr(rot90(indata(:,:,temp_cntr)));
        %             end
        %         end
        
        indata = orintate_image(indata,flag_orient);
        %{
        switch flag_orient
            case 0
            case 1

            for temp_cntr=1:size(indata,3)
                indata(:,:,temp_cntr)=rot90(indata(:,:,temp_cntr),-1);
            end
            case 2
            for temp_cntr=1:size(indata,3)
                indata(:,:,temp_cntr)=fliplr(rot90(indata(:,:,temp_cntr)));
            end
            case 3
            for temp_cntr=1:size(indata,3)
                indata(:,:,temp_cntr)=rot90(indata(:,:,temp_cntr),2);
            end
            case 4
            for temp_cntr=1:size(indata,3)
                indata(:,:,temp_cntr)=flipud(rot90(indata(:,:,temp_cntr)));
            end
            case 5
            for temp_cntr=1:size(indata,3)
                indata(:,:,temp_cntr)=fliplr((indata(:,:,temp_cntr)));
            end
            case 6
            for temp_cntr=1:size(indata,3)
                indata(:,:,temp_cntr)=flipud((indata(:,:,temp_cntr)));
            end
            case 7
            for temp_cntr=1:size(indata,3)
                indata(:,:,temp_cntr)=fliplr(rot90(indata(:,:,temp_cntr),-1));
            end
            case 8
            for temp_cntr=1:size(indata,3)
                indata(:,:,temp_cntr)=(rot90(indata(:,:,temp_cntr)));
            end
        end
        %}

    catch
        disp('Problem with reading dicominfo from '), d
        keyboard
    end
    
    n = n+1;
    if n == 1 || isempty(find(slice_loc == d.SliceLocation))
        slice_loc = [slice_loc,d.SliceLocation];
    else
        break;
    end;
end;
%{
slice_loc2use = slice_loc(1);   % this is like 20.582

slice_match_indices = [];
disp('Please wait... Sorting through the template file directory')
%h1 = waitbar(0,'Please wait... Sorting through the template file directory');
clear sortingValue

for i = 1:length(h)
    %    waitbar(i/length(h),h1);
    d = dicominfo(strcat(templatedirtxt,'/',h(i).name));
%     if d.SliceLocation == slice_loc2use
        slice_match_indices = [slice_match_indices,i];
        sortingValue(i) = str2num(d.AcquisitionTime);
        %sortingValue(i) = d.InstanceNumber;
%     else
%         disp('Note that a slice location does not match!! Likely an incorrect dicom file in this folder ')
%         disp('Need to resolve this')
%         keyboard;
%     end;
end;
%close(h1);
[sortingValue IX] = sort(sortingValue);
h = h(IX);
h = h(2:end);
slice_match_indices = slice_match_indices(IX);
slice_match_indices = slice_match_indices(2:end);

% -- Start the dicom write process -- %
if length(slice_match_indices) ~= size(indata,3)
    %errordlg('Data size (3rd dimension) mismatch, please check input dataset! Will continue','Dataset size difference');
    disp('Data size (3rd dimension) mismatch, please check input dataset! Will continue ');
end;
%}
%h2 = waitbar(0,'Please wait... Writing out the new dicom files');
disp('Please wait... Writing out the new dicom files')

%%%% to find max intensity in the sequence
%max_int=0;
%for k = 1:size(indata,3)
%    templateFilename = strcat(templatedirtxt,'/',h(k).name);
%    tmpg=dicominfo(templateFilename);
%
%if(double(tmpg.LargestImagePixelValue)>max_int)
%        max_int=double(tmpg.LargestImagePixelValue);
%    end
%end;
%indata=indata*(max_int/max(indata(:)));


for k = 1:size(indata,3)
    %if mod(k,10) == 0
    %    fprintf([num2str(round(k/size(indata,3)*10000)/100),'%%...\n'])
    %end
    templateFilename = strcat(templatedirtxt,'/',h(1).name);
    
    addDicomHeader(templateFilename, indata(:,:,k), seriesnum, ...
        seriesdesc, k, outname,cardiac_phase);
    % waitbar(k/(size(indata,3)), h2);
end;
%close(h2);
return;










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
function addDicomHeader(templateFilename, data_in, seriesNum, seriesDesc, instanceNum, outname,cardiac_phase)

%   Copyright UCAIR, University of Utah, 2007.
%   Author: Edward DiBella
%   Last Modified by: Sathya Vijayakumar Dec 6th 2007.


md=dicominfo(templateFilename);

md.SeriesDescription = seriesDesc;

if(strcmp(cardiac_phase,'diastole')==1)
md.ProtocolName = '3D_perfusion_diastole';
elseif(strcmp(cardiac_phase,'systole')==1)
md.ProtocolName = '3D_perfusion_systole';
end

    
md.SeriesNumber = seriesNum;

data=data_in;
flatdata=reshape(data,1,prod(size(data)));


md.Height=size(data,1);  md.Width=size(data,2);
outdata=double(data);

range=[];
if (isempty(range))
    range = [min(flatdata) max(flatdata)];
    scale = min(65535,max(flatdata));
end;
badData = find(~isfinite(outdata));

%scale= 65535/max(flatdata);
% scale= 511/max(flatdata);  % more like Siemens range, trying this 7/15/09 EVRD to see if quantization noise
%scale= max_int/max(flatdata);
scale=1;
outdata = uint16(scale*(outdata));

outdata(badData) = range(1);

md.BitDepth = 16;
md.BitsAllocated = 16;
md.BitsStored = 16;

md.SmallestImagePixelValue = scale*min(flatdata);
md.LargestImagePixelValue = scale*max(flatdata);
md.InstanceNumber = instanceNum;


dicomwrite(squeeze(outdata),strcat(outname, num2strWithZeros(instanceNum),'.dcm') ,md);    % was outdata'



function outnum=num2strWithZeros(num)
if (num<=9)
    outnum=strcat('000',num2str(num));
elseif ((9<num) && num<=99)
    outnum=strcat('00',num2str(num));
elseif ((99<num) && num<=999)
    outnum=strcat('0',num2str(num));
elseif ((999<num) && num<=9999)
    outnum=num2str(num);
else
    disp('Problem with out of range number!!!')
end
return



function [MIDinfolookup,SeriesNumberinfolookup] = findMIDrelations()

myprocessingDir = pwd;


try
    MIDinfolookup = containers.Map('Keytype','int32','ValueType','any');
    SeriesNumberinfolookup = containers.Map('Keytype','int32','ValueType','any');
catch
    clear temp
    temp.bla = 0;
    MIDinfolookup = containers.Map({int32(1)},{temp}); remove(MIDinfolookup,int32(1));
    SeriesNumberinfolookup = containers.Map({int32(1)},{temp}); remove(SeriesNumberinfolookup,int32(1));
end
clc
disp('-----------------------');
if(~strfind(myprocessingDir,'mat_files'))
    disp('I must have the dicom folder present, as well as the RawData folder')
    return;
end

tmp=pwd;
tmp2=strfind(tmp,'ReconData');
raid1path=tmp(1:(tmp2-1));
cd(raid1path)


if(~exist('RawData'))
    cd(myprocessingDir);
    disp('I cannot continue without the RawData folder present.  Please link to it');
    return
end
if(~exist('DicomData'))
    cd(myprocessingDir);
    disp('I cannot continue without the DicomData folder present.  Please link to it or create it');
    return;
else
    cd('DicomData');
    dicomDataFolder = pwd;
    cd(raid1path)
    
end
cd('RawData');
rawDataFolder = pwd;
rawFiles = dir('*MID*.mat');

%filter out all dat files that don't have a mat files that was created, or
%in other words, kill all dat files that haven't been converted to mat
%files
CollectedFiles = [];
OnesToKill = [];
for i=1:length(rawFiles)
    A = sscanf(rawFiles(i).name,'meas_MID%d_%s');
    if(isempty(A)), continue, end
    MID = A(1);
    sisterFiles = dir(['*MID' num2str(MID) '*']);
    hasSister = 0;
    for j=1:length(sisterFiles)
        if(~strcmp(sisterFiles(j).name,rawFiles(i).name))
            hasSister = 1;
        end
    end
    if(~hasSister)
        OnesToKill = [OnesToKill i];
    end
end

if(length(OnesToKill) == length(rawFiles))
    disp('I don''t think you''ve run any conversion script to convert the dat files to mat files');
    disp('I''ll assume that all the first MID number with it''s corresponding SeriesDescription');
    disp('gets to claim the dicoms that have that SeriesDescription');
    disp('And that all subsequent dat files will be thrown away');
else
    disp('There were dat files that didn''t get converted.  They are:');
    for i=1:length(OnesToKill)
        disp(['- ' rawFiles(OnesToKill(i)).name]);
    end
    rawFiles(OnesToKill) = [];
end




disp('-----------------------');
disp('MID: SeriesDescription');
disp('-Series Numbers found');
disp(' ');
SeriesNumbersTaken = zeros(30,1);
for i=length(rawFiles):-1:1
    cd(myprocessingDir);
    A = sscanf(rawFiles(i).name,'meas_MID%d_%s');
    if(isempty(A))
        disp(['Could not parse(' rawFiles(i).name ')']);
        disp('Please construct a struct with the following parameters: ');
        disp('mystruct.SeriesNumber = 1;');
        disp('mystruct.SliceNumber = 1;');
        disp('mystruct.SeriesDescription = ''bla'';');
        A = sscanf(rawFiles(i).name,'meas_MID%d');
        if(~isempty(A))
            disp(['MIDinfolookup(' num2str(A(1)) ') =  mystruct;']);
        else
            disp('MIDinfolookup(someMID) = mystruct;' );
        end
        continue;
    end
    clear mystruct
    MID = A(1);
    
    %extract the seriesDescripotion
    mystruct.SeriesDescription = ['' A(2:end)'];
    FIDindex = strfind(mystruct.SeriesDescription,'_FID');
    mystruct.SeriesDescription = mystruct.SeriesDescription(1:(FIDindex-1));
    description = mystruct.SeriesDescription;
    disp([num2str(MID) ' : ' mystruct.SeriesDescription]);
    
    
    cd(dicomDataFolder);
    searchString = [mystruct.SeriesDescription '(*'];
    searchString = strrep(searchString,'_','*');
    searchString = strrep(searchString,'.','*');
    searchString = strrep(searchString,'/','*');
    candidateDicomFolders = dir(searchString);
    if(length(candidateDicomFolders) == 0)
        disp(['      Cannot find any dicoms with Description: ' mystruct.SeriesDescription ]);
        continue;
    else
        tempg=candidateDicomFolders(1).name;
        description=tempg(1:end-4);
    end
    candidateSeriesNumbers = [];
    for j=1:length(candidateDicomFolders)
        cd(candidateDicomFolders(j).name);
        dicoms = dir('*.dcm');
        dicomHeader = dicominfo(dicoms(1).name);
        cd(dicomDataFolder);
        candidateDicomFolders(j).timestamp = str2num(dicomHeader.AcquisitionTime);
        candidateDicomFolders(j).seriesNumber = dicomHeader.SeriesNumber;
        candidateSeriesNumbers = horzcat(candidateSeriesNumbers,dicomHeader.SeriesNumber);
    end
    [candidateSeriesNumbers,IX] = sort(candidateSeriesNumbers);
    candidateDicomFolders = candidateDicomFolders(IX);
    j = length(candidateDicomFolders);
    collectedSeriesNumbers = [];
    %while
    %(
    %either the next series number that is to be claimed leads to an
    %invalid index in the array that holds which series numbers have been
    %claimed
    %   or
    %that series number has not been claimed
    %)
    %   and
    %we havent' walked off the edge of the list of series numbers to be
    %claimed (this should be the same as the first test)
    
    while(0 < j)
        while(length(candidateSeriesNumbers) >= j && ...
                length(SeriesNumbersTaken) >= candidateSeriesNumbers(j) && ...
                SeriesNumbersTaken(candidateSeriesNumbers(j)) && ...
                j > 0)
            j = j - 1;
            if(j == 0) % GA added Apr 13, 2011
                break;
            end % end
        end
        if(j == 0)
            break;
        end
        collectedSeriesNumbers = [collectedSeriesNumbers candidateSeriesNumbers(j)];
        SeriesNumbersTaken(candidateSeriesNumbers(j)) = 1;
        j = j - 1;
        %if either this is the last one or the the difference to the next
        %timestamp is more than 5 seconds we break due to that being
        %another injection
        if(j == 0 || abs(candidateDicomFolders(j).timestamp - candidateDicomFolders(j+1).timestamp) > 5)
            break;
        end
    end
    
    mystruct.SeriesDescription = description; % GA
    mystruct.SeriesNumbers = sort(collectedSeriesNumbers);
    disp(['      - [' num2str(mystruct.SeriesNumbers) ']']);
    MIDinfolookup(MID) = mystruct;
    clear mystruct
    mystruct.MID = MID;
    mystruct.SeriesDescription = description;
    for slice = 1:length(collectedSeriesNumbers)
        SeriesNumber = collectedSeriesNumbers(slice);
        mystruct.slice = slice;
        if(isKey(SeriesNumberinfolookup,SeriesNumber))
            disp(['Colliding data for series number: ' num2str(SeriesNumber)]);
            disp('Skipping');
            continue;
        end
        SeriesNumberinfolookup(SeriesNumber) = mystruct;
    end
end

cd(myprocessingDir);

save('MIDinfolookup.mat','MIDinfolookup','SeriesNumberinfolookup');

return;

