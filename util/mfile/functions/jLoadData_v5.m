% --------------------------------------------------------------------------------
% [Data,Theta]=jLoadData(FileInfo)
% Jason Mendes - UCAIR - University of Utah - 2015
% Updated July 2016 for multicontrast AIF
% --------------------------------------------------------------------------------
% Load data from measurement file.
% --------------------------------------------------------------------------------
function [Data]=jLoadData_v5(FileInfo)

% If FileInfo is filename then open file
if (~isstruct(FileInfo))
    FileInfo=jGetFileInfo(FileInfo);
    %FileInfo=jGetFileInfo('/v/raid1a/gadluru/MRIdata/Cardiac/Prisma/P033117/RawData/meas_MID00073_FID36826_RadialDCE_PerfusionStress.dat');
    CloseFileOnExit=1;
else
    CloseFileOnExit=0;
end

% Determine configuration and initialize
SatNormal=0;
SatEfficient=1;
Data=[];
N=length(FileInfo.Protocol);
Data.Info.NumberOfSets3D=FileInfo.Protocol(N).sWipMemBlock.alFree(2);
Data.Info.SaturationRecoveryTime3D=FileInfo.Protocol(N).sWipMemBlock.alFree(3);
Data.Info.SetSpacing3D=FileInfo.Protocol(N).sWipMemBlock.alFree(4);
Data.Info.NumberOfLinesAIF=FileInfo.Protocol(N).sWipMemBlock.alFree(5);
Data.Info.SliceThicknessAIF=FileInfo.Protocol(N).sWipMemBlock.alFree(6);
Data.Info.SaturationRecoveryTimeAIF=FileInfo.Protocol(N).sWipMemBlock.alFree(7);
Data.Info.NumberOfContrastsAIF=FileInfo.Protocol(N).sWipMemBlock.alFree(8);
Data.Info.NumberOfProtonDensityScans=FileInfo.Protocol(N).sWipMemBlock.alFree(10);
Data.Info.ProtonDensityFlipAngle=FileInfo.Protocol(N).sWipMemBlock.alFree(11);
Data.Info.VariableDensityRatio=FileInfo.Protocol(N).sWipMemBlock.alFree(12)/100;
Data.Info.AsymmetricRatio=FileInfo.Protocol(N).sWipMemBlock.alFree(13)/100;
if (length(FileInfo.Protocol(N).sWipMemBlock.alFree)>13)
    SatType=FileInfo.Protocol(N).sWipMemBlock.alFree(14);
else
    SatType=0;
end
switch SatType
    case SatNormal
        Data.Info.SequenceType='Normal: AIF before 3D and during 3D SRT';
        Data.Info.NumberOfSetsAIF=Data.Info.NumberOfSets3D+1;
    case SatEfficient
        Data.Info.SequenceType='Compressed: AIF during 3D SRT';
        Data.Info.NumberOfSetsAIF=Data.Info.NumberOfSets3D;
    otherwise
        Data.Info.SequenceType='Unknown';
        Data.Info.NumberOfSetsAIF=1;
end
Data.Info.RepetitionTime3D=FileInfo.Protocol(N).alTR(2)/1e3;
Data.Info.RepetitionTimeAIF=FileInfo.Protocol(N).alTR(3)/1e3;
Data.Info.AcquireTime=FileInfo.Protocol(N).alTR(1)/1e3;
Data.Info.EchoTime3D=FileInfo.Protocol(N).alTE(1)/1e3;
Data.Info.EchoTimeAIF=FileInfo.Protocol(N).alTE(2)/1e3;
Data.Info.EchoSpacingAIF=FileInfo.Protocol(N).alTE(3)/1e3;
try
    Data.Info.TriggerDelay=FileInfo.Protocol(N).alTD(1)/1e3;
catch
    Data.Info.TriggerDelay=0;
end
Data.Info.FlipAngle=FileInfo.Protocol(N).adFlipAngleDegree;
Data.Info.Partitions3D=FileInfo.Protocol(N).sKSpace.lPartitions;
Data.Info.LinesInCentralPartition=FileInfo.Protocol(N).sKSpace.lRadialViews;
Data.Info.ReadoutSamples=2*FileInfo.Protocol(N).sKSpace.lBaseResolution;
Data.Info.NumberOfChannels=FileInfo.Protocol(N).lNumberOfChannels;
Data.Info.NumberOfMeasurements=FileInfo.Protocol(N).lRepetitions+1;
for BlockIndex=1:length(FileInfo.Protocol(N).OffsetImage)
    Header=jGetHeader(FileInfo,BlockIndex);
    if (Header.Repetition>0)
        break;
    end
end
if (BlockIndex<length(FileInfo.Protocol(N).OffsetImage))
    BlocksPerMeasurement=BlockIndex-1;
    Data.Info.NumberOfMeasurements=floor(length(FileInfo.Protocol(N).OffsetImage)/BlocksPerMeasurement);
else
    error('jLoadData_v5(): Unable to find end of first data block');
end
if (BlocksPerMeasurement~=(length(FileInfo.Protocol(N).OffsetImage)/Data.Info.NumberOfMeasurements))
    fprintf('Warning: Data file might be corrupt\n')
    fprintf('Loading in %d of %d measurements\n',Data.Info.NumberOfMeasurements,FileInfo.Protocol(N).lRepetitions+1);
end
Data.Info.NumberOfLines3D=(BlocksPerMeasurement-Data.Info.NumberOfLinesAIF*Data.Info.NumberOfContrastsAIF*Data.Info.NumberOfSetsAIF)/Data.Info.NumberOfSets3D;

% Determine data locations
BlocksPerAIF=Data.Info.NumberOfLinesAIF*Data.Info.NumberOfContrastsAIF;
if (SatType==SatNormal)
    BlockOffsetAIF=[0 BlocksPerAIF 2*BlocksPerAIF+Data.Info.NumberOfLines3D];
    BlockOffset3D=[2*BlocksPerAIF 3*BlocksPerAIF+Data.Info.NumberOfLines3D];
    BlockOffsetAIF=BlockOffsetAIF(1:(Data.Info.NumberOfSets3D+1));
    BlockOffset3D=BlockOffset3D(1:(Data.Info.NumberOfSets3D));
    BlockNameAIF(1).String='Normal';
    BlockNameAIF(2).String='Systole';
    BlockNameAIF(3).String='Diastole';
    BlockNameAIF(4).String='Unknown';
else
    BlockOffsetAIF=[0 BlocksPerAIF+Data.Info.NumberOfLines3D];
    BlockOffset3D=[BlocksPerAIF 2*BlocksPerAIF+Data.Info.NumberOfLines3D];
    BlockOffsetAIF=BlockOffsetAIF(1:Data.Info.NumberOfSets3D);
    BlockOffset3D=BlockOffset3D(1:Data.Info.NumberOfSets3D);
    BlockNameAIF(1).String='Systole';
    BlockNameAIF(2).String='Diastole';
    BlockNameAIF(3).String='Unknown';
    BlockNameAIF(4).String='Unknown';
end
BlockName3D(1).String='Systole';
BlockName3D(2).String='Diastole';
BlockName3D(3).String='Unknown';
BlockOffset3D=BlockOffset3D(1:Data.Info.NumberOfSets3D);
% Read in AIF data
for SetIndex=1:length(BlockOffsetAIF)
    CurrentData=single(complex(0,0));
    CurrentData(Data.Info.ReadoutSamples,Data.Info.NumberOfLinesAIF,1,Data.Info.NumberOfChannels,Data.Info.NumberOfMeasurements,Data.Info.NumberOfContrastsAIF)=single(complex(0,0));
    CurrentData(:)=single(complex(0,0));
    Theta=zeros(Data.Info.NumberOfMeasurements,Data.Info.NumberOfLinesAIF);
    fprintf('Reading in %s AIF...',BlockNameAIF(min(SetIndex,length(BlockNameAIF))).String);
    for MeasIndex=1:Data.Info.NumberOfMeasurements
        for BlockIndex=1:BlocksPerAIF
            [ComplexData,Header]=jGetDataAndHeader(FileInfo,BlockOffsetAIF(SetIndex)+BlockIndex+(MeasIndex-1)*BlocksPerMeasurement);
            Ky=Header.Line+1;
            Kx=(1:Header.Samples)-Header.Samples+Data.Info.ReadoutSamples;
            Theta(MeasIndex,Ky)=Header.IceParam(5)/100;
            Header.Echo=rem(BlockIndex-1,Data.Info.NumberOfContrastsAIF);
            if (rem(Header.Echo,2))
                CurrentData(Kx,Ky,1,:,MeasIndex,Header.Echo+1)=ComplexData(Header.Samples:-1:1,:,:,:);
            else
                CurrentData(Kx,Ky,1,:,MeasIndex,Header.Echo+1)=ComplexData;
            end
        end
    end
    eval(sprintf('Data.AIF_%s=permute(CurrentData,[1 2 3 5 6 4]);',BlockNameAIF(min(SetIndex,length(BlockNameAIF))).String));
    eval(sprintf('Data.AIF_%s_Theta=Theta;',BlockNameAIF(min(SetIndex,length(BlockNameAIF))).String));
    fprintf('Done\n');
end

% Read in 3D data
Data.Info.TriggerTime=zeros(1,Data.Info.NumberOfMeasurements);
for SetIndex=1:length(BlockOffset3D)
    CurrentData=single(complex(0,0));
    CurrentData(Data.Info.ReadoutSamples,Data.Info.LinesInCentralPartition,Data.Info.Partitions3D,Data.Info.NumberOfChannels,Data.Info.NumberOfMeasurements)=single(complex(0,0));
    CurrentData(:)=single(complex(0,0));
    Theta=zeros(Data.Info.NumberOfMeasurements,Data.Info.LinesInCentralPartition,Data.Info.Partitions3D);
    fprintf('Reading in %s 3D...',BlockName3D(min(SetIndex,length(BlockName3D))).String);
    for MeasIndex=1:Data.Info.NumberOfMeasurements
        for BlockIndex=1:Data.Info.NumberOfLines3D
            [ComplexData,Header]=jGetDataAndHeader(FileInfo,BlockOffset3D(SetIndex)+BlockIndex+(MeasIndex-1)*BlocksPerMeasurement);
            Ky=Header.Line+1;
            Kz=Header.Partition-Header.KspacePartitionCenter+floor(Data.Info.Partitions3D/2)+1;
            Kx=(1:Header.Samples)-Header.Samples+Data.Info.ReadoutSamples;
            Theta(MeasIndex,Ky,Kz)=Header.IceParam(5)/100;
            CurrentData(Kx,Ky,Kz,:,MeasIndex)=ComplexData;
            if (BlockIndex==1)
                Data.Info.TriggerTime(MeasIndex)=Header.Time;
            end
        end
    end
    eval(sprintf('Data.%s=permute(CurrentData,[1 2 3 5 4]);',BlockName3D(min(SetIndex,length(BlockName3D))).String));
    eval(sprintf('Data.%s_Theta=Theta;',BlockName3D(min(SetIndex,length(BlockName3D))).String));
    fprintf('Done\n');
end
Data.Info.TriggerTime=Data.Info.TriggerTime-Data.Info.TriggerTime(1);

% Close file if opned in this function
if (CloseFileOnExit)
    fclose(FileInfo.FileId);
end


% --------------------------------------------------------------------------------
% [ComplexData,Header]=jGetDataAndHeader(FileInfo,BlockIndex,SequenceIndex)
% Jason Mendes - UCAIR - University of Utah - 2015
% --------------------------------------------------------------------------------
% This function reads in the data and header from a Siemens .dat raw measurement
% data file.
% --------------------------------------------------------------------------------
function [ComplexData,Header]=jGetDataAndHeader(FileInfo,BlockIndex,SequenceIndex)

if (exist('BlockIndex','var')~=1)
   BlockIndex=1;
end
if (exist('SequenceIndex','var')~=1)
   SequenceIndex=length(FileInfo.Header);
end
Header=jGetHeader(FileInfo,BlockIndex,SequenceIndex);
ComplexData=jGetData(FileInfo,BlockIndex,SequenceIndex);


% --------------------------------------------------------------------------------
% ComplexData=jGetData(FileInfo,BlockIndex,SequenceIndex)
% Jason Mendes - UCAIR - University of Utah - 2015
% --------------------------------------------------------------------------------
% This function reads in the data from a Siemens .dat raw measurement data file.
% --------------------------------------------------------------------------------
function ComplexData=jGetData(FileInfo,BlockIndex,SequenceIndex)

if (exist('BlockIndex','var')~=1)
   BlockIndex=1;
end
if (exist('SequenceIndex','var')~=1)
   SequenceIndex=length(FileInfo.Header);
end

% VB style
if (FileInfo.IsVB)
   fseek(FileInfo.FileId,FileInfo.Protocol(SequenceIndex).OffsetImage(BlockIndex),'bof');
   Data=fread(FileInfo.FileId,3,'uchar=>uchar');
   LengthDMA=double(typecast([Data(1:3);0],'uint32'));
   fseek(FileInfo.FileId,FileInfo.Protocol(SequenceIndex).OffsetImage(BlockIndex),'bof');
   Data=fread(FileInfo.FileId,LengthDMA,'uchar=>uchar');   
   Nc=FileInfo.Protocol(SequenceIndex).lNumberOfChannels;
   ChannelSize=LengthDMA/Nc;
   Nx=(ChannelSize-128)/4/2;
   ComplexData=single(complex(0,0));
   ComplexData(Nx,1,1,Nc)=0;
   for ChannelIndex=1:Nc
      RawData=typecast(Data((129:ChannelSize)+(ChannelIndex-1)*ChannelSize),'single');
      ComplexData(:,1,1,ChannelIndex)=complex(RawData(1:2:2*Nx),RawData(2:2:2*Nx));
   end
   
% VD style
elseif (FileInfo.IsVD)
   fseek(FileInfo.FileId,FileInfo.Protocol(SequenceIndex).OffsetImage(BlockIndex),'bof');
   Data=fread(FileInfo.FileId,3,'uchar=>uchar');
   LengthDMA=double(typecast([Data(1:3);0],'uint32'));
   fseek(FileInfo.FileId,FileInfo.Protocol(SequenceIndex).OffsetImage(BlockIndex),'bof');
   Data=fread(FileInfo.FileId,LengthDMA,'uchar=>uchar');   
   Nc=FileInfo.Protocol(SequenceIndex).lNumberOfChannels;
   ChannelSize=(LengthDMA-192)/Nc;
   Nx=(ChannelSize-32)/4/2;
   ComplexData=single(1i);
   ComplexData(Nx,1,1,Nc)=0;
   for ChannelIndex=1:Nc
      RawData=typecast(Data((33:ChannelSize)+(ChannelIndex-1)*ChannelSize+192),'single');
      ComplexData(:,1,1,ChannelIndex)=complex(RawData(1:2:2*Nx),RawData(2:2:2*Nx));
   end
   
 % Unknown type
else
   ComplexData=[];
end


% --------------------------------------------------------------------------------
%  Header=jGetHeader(FileInfo,BlockIndex,SequenceIndex)
% Jason Mendes - UCAIR - University of Utah - 2015
% --------------------------------------------------------------------------------
% This function reads in just the header data from a Siemens .dat raw measurement
% data file.
% --------------------------------------------------------------------------------
function Header=jGetHeader(FileInfo,BlockIndex,SequenceIndex)

if (exist('BlockIndex','var')~=1)
   BlockIndex=1;
end
if (exist('SequenceIndex','var')~=1)
   SequenceIndex=length(FileInfo.Header);
end

% VB style
if (FileInfo.IsVB)
   fseek(FileInfo.FileId,FileInfo.Protocol(SequenceIndex).OffsetImage(BlockIndex),'bof');
   Data=fread(FileInfo.FileId,128,'uchar=>uchar');
   Header.LengthDMA=double(typecast([Data(1:3);0],'uint32'));
   Header.Id=double(typecast(Data(5:8),'int32'));
   Header.Scan=double(typecast(Data(9:12),'uint32'));
   Header.Time=0.0025*double(typecast(Data(13:16),'uint32'));
   Header.Pmu=0.0025*double(typecast(Data(17:20),'uint32'));
   Header.Mask=double(typecast(Data(21:28),'uint32'));
   Header.Samples=double(typecast(Data(29:30),'uint16'));
   Header.Channels=double(typecast(Data(31:32),'uint16'));
   Header.Line=double(typecast(Data(33:34),'uint16'));
   Header.Acquisition=double(typecast(Data(35:36),'uint16'));
   Header.Slice=double(typecast(Data(37:38),'uint16'));
   Header.Partition=double(typecast(Data(39:40),'uint16'));
   Header.Echo=double(typecast(Data(41:42),'uint16'));
   Header.Phase=double(typecast(Data(43:44),'uint16'));
   Header.Repetition=double(typecast(Data(45:46),'uint16'));
   Header.Set=double(typecast(Data(47:48),'uint16'));
   Header.Seg=double(typecast(Data(49:50),'uint16'));
   Header.Ida=double(typecast(Data(51:52),'uint16'));  
   Header.Idb=double(typecast(Data(53:54),'uint16')); 
   Header.Idc=double(typecast(Data(55:56),'uint16')); 
   Header.Idd=double(typecast(Data(57:58),'uint16')); 
   Header.Ide=double(typecast(Data(59:60),'uint16')); 
   Header.PreCutoff=double(typecast(Data(61:62),'uint16'));
   Header.PostCutoff=double(typecast(Data(63:64),'uint16'));
   Header.KspaceReadCenter=double(typecast(Data(65:66),'uint16'));
   Header.CoilSelect=double(typecast(Data(67:68),'uint16'));
   Header.ReadOffset=double(typecast(Data(69:72),'single'));
   Header.TimeSinceLastRF=0.0025*double(typecast(Data(73:76),'uint32'));
   Header.KspaceLineCenter=double(typecast(Data(77:78),'uint16'));
   Header.KspacePartitionCenter=double(typecast(Data(79:80),'uint16'));
   Header.IceParam=double(typecast(Data(81:96),'uint16'));
   Header.Sag=double(typecast(Data(97:100),'single'));
   Header.Cor=double(typecast(Data(101:104),'single'));
   Header.Tra=double(typecast(Data(105:108),'single'));
   Header.Rotate=double(typecast(Data(109:124),'single'));
   Header.Channel=double(typecast(Data(125:126),'uint16'));
   Header.PTabPos=double(typecast(Data(127:128),'uint16'));
   Header.IsPMU=min(bitand(Header.Mask(1),2^5),1);
   Header.Reflect=min(bitand(Header.Mask(1),2^24),1);
   Header.PhaseCorrect=min(bitand(Header.Mask(1),2^21),1);
   Header.PhaseFFT=min(bitand(Header.Mask(1),2^18),1);
   Header.PartitionFFT=min(bitand(Header.Mask(1),2^16),1);
   Header.LastScanInMeasurment=min(bitand(Header.Mask(1),2^11),1);
   Header.GainSwitch=min(bitand(Header.Mask(1),2^10),1);
   Header.CalibrateKWIC=min(bitand(Header.Mask(2),2^3),1);
   
% VD style
elseif (FileInfo.IsVD)
   fseek(FileInfo.FileId,FileInfo.Protocol(SequenceIndex).OffsetImage(BlockIndex),'bof');
   Data=fread(FileInfo.FileId,192,'uchar=>uchar');
   Header.LengthDMA=double(typecast([Data(1:3);0],'uint32'));
   Header.Id=double(typecast(Data(5:8),'int32'));
   Header.Scan=double(typecast(Data(9:12),'uint32'));
   Header.Time=0.0025*double(typecast(Data(13:16),'uint32'));
   Header.Pmu=0.0025*double(typecast(Data(17:20),'uint32'));
   Header.System=double(typecast(Data(21:22),'uint16'));
   Header.PTabDelay=double(typecast(Data(23:24),'uint16'));
   Header.PTabPosX=double(typecast(Data(25:28),'int32'));
   Header.PTabPosY=double(typecast(Data(29:32),'int32'));
   Header.PTabPosZ=double(typecast(Data(33:36),'int32'));
   Header.Reserve=double(typecast(Data(37:40),'uint32'));
   Header.Mask=double(typecast(Data(41:48),'uint32'));
   Header.Samples=double(typecast(Data(49:50),'uint16'));
   Header.Channels=double(typecast(Data(51:52),'uint16'));
   Header.Line=double(typecast(Data(53:54),'uint16'));
   Header.Acquisition=double(typecast(Data(55:56),'uint16'));
   Header.Slice=double(typecast(Data(57:58),'uint16'));
   Header.Partition=double(typecast(Data(59:60),'uint16'));
   Header.Echo=double(typecast(Data(61:62),'uint16'));
   Header.Phase=double(typecast(Data(63:64),'uint16'));
   Header.Repetition=double(typecast(Data(65:66),'uint16'));
   Header.Set=double(typecast(Data(67:68),'uint16'));
   Header.Seg=double(typecast(Data(69:70),'uint16'));
   Header.Ida=double(typecast(Data(71:72),'uint16'));
   Header.Idb=double(typecast(Data(73:74),'uint16'));
   Header.Idc=double(typecast(Data(75:76),'uint16'));
   Header.Idd=double(typecast(Data(77:78),'uint16'));
   Header.Ide=double(typecast(Data(79:80),'uint16'));
   Header.PreCutoff=double(typecast(Data(81:82),'uint16'));
   Header.PostCutoff=double(typecast(Data(83:84),'uint16'));
   Header.KspaceReadCenter=double(typecast(Data(85:86),'uint16'));
   Header.CoilSelect=double(typecast(Data(87:88),'uint16'));
   Header.ReadOffset=double(typecast(Data(89:92),'single'));
   Header.TimeSinceLastRF=0.0025*double(typecast(Data(93:96),'uint32'));
   Header.KspaceLineCenter=double(typecast(Data(97:98),'uint16'));
   Header.KspacePartitionCenter=double(typecast(Data(99:100),'uint16'));
   Header.Sag=double(typecast(Data(101:104),'single'));
   Header.Cor=double(typecast(Data(105:108),'single'));
   Header.Tra=double(typecast(Data(109:112),'single'));
   Header.Rotate=double(typecast(Data(113:128),'single'));
   Header.IceParam=double(typecast(Data(129:176),'uint16'));
   Header.ReservePara=double(typecast(Data(177:184),'uint16'));
   Header.ApplicationCounter=double(typecast(Data(185:186),'uint16'));
   Header.ApplicationMask=double(typecast(Data(187:188),'uint16'));
   Header.CRC=double(typecast(Data(189:192),'uint32'));
   Header.IsPMU=min(bitand(Header.Mask(1),2^5),1);
   Header.Reflect=min(bitand(Header.Mask(1),2^24),1);
   Header.PhaseCorrect=min(bitand(Header.Mask(1),2^21),1);
   Header.PhaseFFT=min(bitand(Header.Mask(1),2^18),1);
   Header.PartitionFFT=min(bitand(Header.Mask(1),2^16),1);
   Header.LastScanInMeasurement=min(bitand(Header.Mask(1),2^11),1);
   Header.GainSwitch=min(bitand(Header.Mask(1),2^10),1);
   Header.CalibrateKWIC=min(bitand(Header.Mask(2),2^3),1);

 % Unknown type
else
   Header=[];
end


% --------------------------------------------------------------------------------
% FileInfo=jGetFileInfo(FileName)
% Jason Mendes - UCAIR - University of Utah - 2015
% --------------------------------------------------------------------------------
% This function opens the .dat file and retrieves information about the sequence
% data stored.
% --------------------------------------------------------------------------------
function FileInfo=jGetFileInfo(FileName)

% Prompt for filename if one is not specified
if ((exist('FileName')~=1)||isempty(FileName))
	[file_name,path_name]=uigetfile('*.dat','Select Raw Data File');
	FileName=sprintf('%s%s',path_name,file_name);
end

% Open the .dat file
FileInfo=[];
FileInfo.FileId=fopen(FileName,'r');
if (FileInfo.FileId<0)
   error('jGetFileInfo(): Could not open file specified');
end
fseek(FileInfo.FileId,0,'eof');
FileSize=ftell(FileInfo.FileId);
FileInfo.FileName=FileName;

% Find BaselineString
fseek(FileInfo.FileId,0,'bof');
while (ftell(FileInfo.FileId)<FileSize)
   NextLine=fgets(FileInfo.FileId);
   if (~isempty(findstr(NextLine,'<ParamString."tMeasuredBaselineString">')))
      Data=jReadSimpleString(FileInfo,NextLine);
      if (~isempty(findstr(lower(Data.String),'n4_v')))
         FileInfo.Version=Data.String;
         break;
      end
   elseif (~isempty(findstr(NextLine,'<ParamString."tBaselineString">')))
      Data=jReadSimpleString(FileInfo,NextLine);
      if (~isempty(findstr(lower(Data.String),'n4_v')))
         FileInfo.Version=Data.String;
         break;
      end
   end
end

% If no version is found try and guess one
if (~isfield(FileInfo,'Version'))
   fseek(FileInfo.FileId,0,'bof');
   Data=fread(FileInfo.FileId,2,'uint32');
   if ((Data(1)<10000)&&(Data(2)<=64))
      FileInfo.Version='BestMatchIsVD';
   else
      FileInfo.Version='BestMatchIsVB';
   end
end
if (~isempty(findstr(lower(FileInfo.Version),'vb')))
   FileInfo.IsVB=1;
   FileInfo.IsVD=0;
elseif (~isempty(findstr(lower(FileInfo.Version),'vd')))
   FileInfo.IsVB=0;
   FileInfo.IsVD=1;
elseif (~isempty(findstr(lower(FileInfo.Version),'ve')))
   FileInfo.IsVB=0;
   FileInfo.IsVD=1;
else
   error('jGetFileInfo(): Unknown version');
end

% Find how many measurements are in the file
if (FileInfo.IsVB)
   fseek(FileInfo.FileId,0,'bof');
   FileInfo.Header.StartPosition=0;
   FileInfo.Header.EndPosition=FileSize;
   FileInfo.Header.DataPosition=fread(FileInfo.FileId,1,'uint32');
elseif (FileInfo.IsVD)
   fseek(FileInfo.FileId,4,'bof');
   NumberOfMeasurements=fread(FileInfo.FileId,1,'uint32');
   for MeasIndex=1:NumberOfMeasurements
      FileInfo.Header(MeasIndex).MeasId=fread(FileInfo.FileId,1,'uint32');
      FileInfo.Header(MeasIndex).FileId=fread(FileInfo.FileId,1,'uint32');
      FileInfo.Header(MeasIndex).StartPosition=fread(FileInfo.FileId,1,'uint64');
      FileInfo.Header(MeasIndex).EndPosition=FileInfo.Header(MeasIndex).StartPosition+fread(FileInfo.FileId,1,'uint64');
      FileInfo.Header(MeasIndex).PatientName=sprintf('%s',fread(FileInfo.FileId,64,'uchar=>char'));
      FileInfo.Header(MeasIndex).ProtocolName=sprintf('%s',fread(FileInfo.FileId,64,'uchar=>char'));
   end
   for MeasIndex=1:NumberOfMeasurements
      if (fseek(FileInfo.FileId,FileInfo.Header(MeasIndex).StartPosition,'bof'))
         error('jGetFileInfo(): Data set not complete');
      end
      FileInfo.Header(MeasIndex).DataPosition=FileInfo.Header(MeasIndex).StartPosition+fread(FileInfo.FileId,1,'uint32');
   end
else
   error('jGetFileInfo(): Unknown file type');
end

% Load all protocols
FileInfo.Protocol=jGetProtocol(FileInfo);


% --------------------------------------------------------------------------------
% FileInfo=jLoadBasicProtocol(FileName)
% Jason Mendes - UCAIR - University of Utah - 2015
% --------------------------------------------------------------------------------
% This function opens the .dat file and retrieves information about the sequence
% protocol from the ASCCONV section and some additional parameters.
% --------------------------------------------------------------------------------
function Protocol=jGetProtocol(FileInfo)

% Prepare to read protocol
Protocol=[];

% There are some useful parameters not in the ASCCONV section
ParameterTag(1).String='<ParamLong."lNoOfPhaseCorrScans">';
ParameterTag(2).String='<ParamLong."relSliceNumber">';
ParameterTag(3).String='<ParamLong."NSetMeas">';
ParameterTag(4).String='<ParamLong."RampupTime">';
ParameterTag(5).String='<ParamLong."FlattopTime">';
ParameterTag(6).String='<ParamLong."DelaySamplesTime">';
ParameterTag(7).String='<ParamLong."RegridMode">';
ParameterTag(8).String='<ParamDouble."ADCDuration">';
ParameterTag(9).String='<ParamLong."alRegridRampupTime">';
ParameterTag(10).String='<ParamLong."alRegridFlattopTime">';
ParameterTag(11).String='<ParamLong."alRegridDelaySamplesTime">';
ParameterTag(12).String='<ParamLong."alRegridMode">';
ParameterTag(13).String='<ParamDouble."aflRegridADCDuration">';
ParameterTag(14).String='<ParamLong."NoOfFourierPartitions">';
ParameterTag(15).String='<ParamString."PatientPosition">';

% Find protocol for each requested sequence
for SearchIndex=1:length(FileInfo.Header)
   % Read in the header
   HeaderStart=FileInfo.Header(SearchIndex).StartPosition;
   HeaderEnd=FileInfo.Header(SearchIndex).DataPosition;
   if (fseek(FileInfo.FileId,HeaderStart,'bof'))
      error('jGetProtocol(): Data header not complete');
   end
   Header=fread(FileInfo.FileId,HeaderEnd-HeaderStart,'char=>char')';
   % Find ASCCONV section
   StartIndex=max(findstr(Header,'ASCCONV BEGIN'));
   EndIndex=max(findstr(Header,'ASCCONV END'));
   if ((~isempty(StartIndex))&&(~isempty(EndIndex)))
      String=Header(StartIndex:EndIndex);
      String(find(String==';'))=':'; % Remove original ';'
      String=regexprep(String,'\n',';\n');
      String=String((find(String==10,1)+1):find(String==10,1,'last'));
      String(find(String<=32))=[]; % Remove spaces
      String(find(String=='_'))=[]; % '_' cause problems
      String=regexprep(String,'[','('); % Matlab uses round brackets for indices
      String=regexprep(String,']',')'); % Matlab uses round brackets for indices
      String=regexprep(String,'"',''''); % Matlab uses round brackets for indices
      String=regexprep(String,'(','(1+'); % Matlab is not zero based
      String=regexprep(String,'0x(\w*)','hex2dec(''$1'')'); % Convert hex values
      LineEnd=[0 find(String==';')];
      for Index=2:length(LineEnd)
         EvalString=String((LineEnd(Index-1)+1):LineEnd(Index));
         if isempty(findstr(lower(EvalString),'attribute'))
            try
               eval(sprintf('Protocol(%d).%s',SearchIndex,EvalString));
            catch
               fprintf('Error adding %s\n',EvalString');
            end
         end
      end
   end   
   % Read in extra parameters
   for Index=1:length(ParameterTag);
      StartIndex=min(findstr(Header,ParameterTag(Index).String));
      EndIndex=find(Header(StartIndex:end)=='}',1);
      if ((~isempty(StartIndex))&&(~isempty(EndIndex)))
         CurrentString=Header((1:EndIndex)+StartIndex);
         NameIndex=find(CurrentString=='"',2);
         ValueIndex=(find(CurrentString=='{',1,'last')+1):(find(CurrentString=='}',1,'last')-1);
         if (length(NameIndex==2)&&length(ValueIndex>0))
            NameString=CurrentString((NameIndex(1)+1):(NameIndex(2)-1));
            ValueString=CurrentString(ValueIndex);
            if (~isempty(findstr(lower(CurrentString),'double')))
               ValueString=ValueString((min(findstr(lower(ValueString),'<precision>'))+12):end);
               if (~isempty(find(ValueString==10,1)))
                   ValueString=ValueString((find(ValueString==10,1)+1):end);
               end
               ValueString=sprintf('%s',ValueString(find(ValueString<=32,1):end));
            elseif (~isempty(findstr(lower(CurrentString),'long')))
               ValueString=sprintf('%s',ValueString);
            elseif (~isempty(findstr(lower(CurrentString),'string')))
               if (length(find(ValueString=='"'))<2)
                  ValueString='''''';
               else
                  ValueString(find(ValueString=='"'))='''';
               end
            else
               % Not supported at this time
               ValueString='[]';
            end
            EvalString=sprintf('Protocol(%d).%s=[%s];',SearchIndex,NameString,ValueString);
            EvalString(find(EvalString==10))=32;
            try
               eval(EvalString);
            catch
               fprintf('Error %s\n',EvalString');
               pause;
            end
         end
      end
   end
   % Compatibility issues
    if (isfield(Protocol(SearchIndex),'ADCDuration'))
        if (~isempty(Protocol(SearchIndex).ADCDuration))
            Protocol(SearchIndex).aflRegridADCDuration=Protocol(SearchIndex).ADCDuration;
        end
    end
    if (isfield(Protocol(SearchIndex),'RampupTime'))
        if (~isempty(Protocol(SearchIndex).RampupTime))
            Protocol(SearchIndex).alRegridRampupTime=Protocol(SearchIndex).RampupTime;
        end
    end
    if (isfield(Protocol(SearchIndex),'FlattopTime'))
        if (~isempty(Protocol(SearchIndex).FlattopTime))
            Protocol(SearchIndex).alRegridFlattopTime=Protocol(SearchIndex).FlattopTime;
        end
    end
    if (isfield(Protocol(SearchIndex),'DelaySamplesTime'))
        if (~isempty(Protocol(SearchIndex).DelaySamplesTime))
            Protocol(SearchIndex).alRegridDelaySamplesTime=Protocol(SearchIndex).DelaySamplesTime;
        end
    end
    if (isfield(Protocol(SearchIndex),'RegridMode'))
        if (~isempty(Protocol(SearchIndex).RegridMode))
            Protocol(SearchIndex).alRegridMode=Protocol(SearchIndex).RegridMode;
        end
    end
    if (isfield(Protocol(SearchIndex),'sWiPMemBlock'))
        if (~isempty(Protocol(SearchIndex).sWiPMemBlock))
            Protocol(SearchIndex).sWipMemBlock=Protocol(SearchIndex).sWiPMemBlock;
        end
    end
   % Fix problem with partition count for 2D
   if (isfield(Protocol(SearchIndex),'sKSpace'))
      if (Protocol(SearchIndex).sKSpace.ucDimension<3)
         Protocol(SearchIndex).sKSpace.lPartitions=1;
      end
   end
   % Single repetition is not recorded
   if (~isfield(Protocol(SearchIndex),'lRepetitions'))
      Protocol(SearchIndex).lRepetitions=0;
   end
   if (isempty(Protocol(SearchIndex).lRepetitions))
       Protocol(SearchIndex).lRepetitions=0;
   end      
    % Number of coils
    if (isfield(Protocol(SearchIndex),'sCoilSelectMeas'))
        Protocol(SearchIndex).asCoilSelectMeas=Protocol(SearchIndex).sCoilSelectMeas;
    end
    if (isfield(Protocol(SearchIndex),'asCoilSelectMeas'))
        if (FileInfo.IsVB)
            Protocol(SearchIndex).lNumberOfChannels=length(Protocol(SearchIndex).asCoilSelectMeas(1).asList);
        elseif (FileInfo.IsVD)
            Protocol(SearchIndex).lNumberOfChannels=length(Protocol(SearchIndex).sCoilSelectMeas(1).aRxCoilSelectData(1).asList);
        end
    end
   % Find location of PMU and image lines
   if (FileInfo.IsVB)
      fseek(FileInfo.FileId,FileInfo.Header(SearchIndex).DataPosition,'bof');
      Data=fread(FileInfo.FileId,3,'uchar=>uchar');
      LengthDMA=double(typecast([Data(1:3);0],'uint32'));
      Protocol(SearchIndex).OffsetPMU=[];
      Protocol(SearchIndex).OffsetImage=(FileInfo.Header(SearchIndex).DataPosition): LengthDMA:(FileInfo.Header(SearchIndex).EndPosition-LengthDMA);
   elseif (FileInfo.IsVD)
      for LoopIndex=1:2
         CountPMU=0;
         CountImage=0;
         CurrentPosition=FileInfo.Header(SearchIndex).DataPosition;
         while (CurrentPosition<=(FileInfo.Header(SearchIndex).EndPosition-192))
            fseek(FileInfo.FileId,CurrentPosition,'bof');
            Data=fread(FileInfo.FileId,3,'uchar=>uchar');
            LengthDMA=double(typecast([Data(1:3);0],'uint32'));
            fseek(FileInfo.FileId,CurrentPosition+40,'bof');
            Data=fread(FileInfo.FileId,4,'uchar=>uchar');
            Mask=double(typecast(Data,'uint32'));
            if (min(bitand(Mask,2^5),1))
               CountPMU=CountPMU+1;
               if (LoopIndex==2)
                  Protocol(SearchIndex).OffsetPMU(CountPMU)=CurrentPosition;
               end
            else
               if (LengthDMA>(Protocol(SearchIndex).lNumberOfChannels*32+192))
                  CountImage=CountImage+1;
                  if (LoopIndex==2)
                     Protocol(SearchIndex).OffsetImage(CountImage)=CurrentPosition;
                  end
               end
            end
            CurrentPosition=CurrentPosition+LengthDMA;
         end
         if (LoopIndex==1)
            Protocol(SearchIndex).OffsetPMU=zeros(1,CountPMU);
            Protocol(SearchIndex).OffsetImage=zeros(1,CountImage);
         end
      end
   end
   % Get readout length
   
end


% --------------------------------------------------------------------------------
% Data=jReadSimpleString(FileInfo,InitialString)
% Jason Mendes - UCAIR - University of Utah - 2014
% Update: Returns both data string name and value (October 2015)
% --------------------------------------------------------------------------------
% This function reads in a string parameter value from a .dat file.
% Returns Data.Name and Data.String containing the string name and value.
% --------------------------------------------------------------------------------
function Data=jReadSimpleString(FileInfo,InitialString)

% Read in data until first closing bracket is found
FinalString=InitialString;
Data=[];
while (isempty(find(FinalString=='}',1)))
   NextString=fgets(FileInfo.FileId,1000);
   FinalString=sprintf('%s %s',FinalString,NextString);
end

% Get field name
ValueIndex=find(FinalString=='"',2,'first');
if (length(ValueIndex)<2)
   DataName=sprintf('UnknownName%d',round(rand(1)*1000));
else
   DataName=FinalString(ValueIndex(1)+1:ValueIndex(2)-1);
end
DataName(find(DataName<33))=[]; % Ignore line feeds, etc...

% Get the values
ValueIndex=find(FinalString=='{',1,'last')+1:find(FinalString=='}',1,'last')-1;
StringText=FinalString(ValueIndex);
ValueIndex=find(StringText=='"',2,'first');
if (length(ValueIndex)<2)
   DataString='Error';
else
   DataString=StringText(ValueIndex(1)+1:ValueIndex(2)-1);
end

% Return values
Data.Name=DataName;
Data.String=DataString;