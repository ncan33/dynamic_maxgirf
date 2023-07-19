% --------------------------------------------------------------------------------
% [DataIR]=jLoadDataIR(FileName,AlignKspaceCenter)
% Jason Mendes - UCAIR - University of Utah - 2018
% --------------------------------------------------------------------------------
% Load data from UCAIR IR sequence. AlignKspaceCenter is used to correct
% for gradient delays. It can be set to 'CalculateAndApply' (default),
% 'CalculateOnly' or 'None'. Align parameters can be found in Info.AlignPixel.
% --------------------------------------------------------------------------------
function [DataIR]=jLoadDataIR(FileName,AlignKspaceCenter)

% Definitions from sequence (UCAIR_UI.h, UCAIR_Perfusion.h and UCAIR_RadialKernel.h)
MultiBandFactor_ID=1;
PulseDuration_ID=2;
MeasurementsPerInversion_ID=3;

% Open file and initialize
fprintf('jLoadDataIR(): Opening file %s\n',FileName);
FileInfo=jGetFileInfo(FileName);
DataIR=[];
Info=[];
Info.Protocol=FileInfo.Protocol(length(FileInfo.Protocol));

% Get parameters from sequence
fprintf('jLoadDataIR(): Saving protocol information\n');
ParameterWIP=zeros(1,14);
ParameterWIP(1:(length(Info.Protocol.sWipMemBlock.alFree)-1))=Info.Protocol.sWipMemBlock.alFree(2:length(Info.Protocol.sWipMemBlock.alFree));
IsGated=Info.Protocol.sPhysioImaging.lSignal1==2;
if (~IsGated)
    Info.Trigger='Not Gated';
else
    Info.Trigger='Gated';
end
Info.NumberOfSlices=length(Info.Protocol.sSliceArray.asSlice);
Info.MultiBandFactor=max(ParameterWIP(MultiBandFactor_ID),1);
Info.MultiBandSets=ceil(Info.NumberOfSlices/Info.MultiBandFactor);
if (Info.Protocol.sTXSPEC.ucRFPulseType==1)
    Info.PulseType='Fast: Single Lobe Sinc';
elseif (Info.Protocol.sTXSPEC.ucRFPulseType==4)
    Info.PulseType='Low Sar: 2 Lobe Asymmetric Sinc';
else
    Info.PulseType='Normal: 4 Lobe Asymmetric Sinc';
end
if (Info.Protocol.sPrepPulses.ucSatRecovery==4)
    Info.SaturationType='UCAIR 6 RF pulses';
elseif (Info.Protocol.sPrepPulses.ucSatRecovery==8)
    Info.SaturationType='Siemens 3 RF pulses';
else
    Info.SaturationType='None';
end
Info.NumberOfMeasurements=Info.Protocol.lRepetitions+1;
Info.ReadoutSamples=2*Info.Protocol.sKSpace.lBaseResolution;
Info.NumberOfChannels=Info.Protocol.lNumberOfChannels;
Info.FlipAngle=Info.Protocol.adFlipAngleDegree;
Info.FOV=Info.Protocol.sSliceArray(1).asSlice(1).dReadoutFOV;
Info.RadialViews=Info.Protocol.sKSpace.lRadialViews;
TR=zeros(1,7);
TR(1:length(Info.Protocol.alTR))=Info.Protocol.alTR;
TE=zeros(1,6);
TE(1:length(Info.Protocol.alTE))=Info.Protocol.alTE;
try
    Info.TriggerDelay_ms=Info.Protocol.alTD(1)*1e-3;
catch
    Info.TriggerDelay=[];
end
Info.TR_ms=TR(1)*1e-3;
Info.TimePerLine_ms=TR(2)*1e-3;
Info.ImageTime_ms=TR(3)*1e-3;
Info.PulseDurationRF_ms=ParameterWIP(PulseDuration_ID);
Info.SliceThickness_mm=Info.Protocol.sSliceArray.asSlice(1).dThickness;
Info.EchoTime_ms=TE(1)*1e-3;
Info.NumberOfInversionTimes=ParameterWIP(MeasurementsPerInversion_ID);
FirstInversionTime_ms=Info.Protocol.alTI(1)*1e-3;
NumberOfLines=length(Info.Protocol.OffsetImage);
if (NumberOfLines~=(Info.RadialViews*Info.MultiBandSets*Info.NumberOfMeasurements))
    fprintf('jLoadDataIR(): >>>> Warning: Data size is not correct <<<<\n');
end
Info.NumberOfInversionPulses=ceil(Info.NumberOfMeasurements/Info.NumberOfInversionTimes);
DataIR.Info=Info;

% Read IR data
fprintf('jLoadDataIR(): Loading IR data\n');
DataIR.Kspace=zeros(DataIR.Info.ReadoutSamples,DataIR.Info.RadialViews,DataIR.Info.NumberOfInversionTimes,DataIR.Info.NumberOfChannels,DataIR.Info.MultiBandSets,DataIR.Info.NumberOfInversionPulses,'single');
DataIR.Theta=zeros(DataIR.Info.RadialViews,DataIR.Info.NumberOfInversionTimes,DataIR.Info.MultiBandSets,DataIR.Info.NumberOfInversionPulses,'single');
DataIR.PhaseRF=zeros(DataIR.Info.RadialViews,DataIR.Info.NumberOfInversionTimes,DataIR.Info.MultiBandSets,DataIR.Info.NumberOfInversionPulses,'single');
InversionTimes=zeros(DataIR.Info.NumberOfInversionTimes,DataIR.Info.MultiBandSets,DataIR.Info.NumberOfInversionPulses,'single');
for LineIndex=1:NumberOfLines
    [ComplexData,Header]=jGetDataAndHeader(FileInfo,LineIndex);
    Kx=(1:Header.Samples)-Header.Samples+DataIR.Info.ReadoutSamples;
    Ky=Header.Line+1;
    InversionTimeIndex=rem(Header.Repetition,DataIR.Info.NumberOfInversionTimes)+1;
    MultiBandSetIndex=Header.IceParam(2)+1;
    InversionPulseIndex=floor(Header.Repetition/DataIR.Info.NumberOfInversionTimes)+1;
    DataIR.Theta(Ky,InversionTimeIndex,MultiBandSetIndex,InversionPulseIndex)=Header.IceParam(5)/100;
    DataIR.PhaseRF(Ky,InversionTimeIndex,MultiBandSetIndex,InversionPulseIndex)=Header.IceParam(7);
    if (Ky==1)
        fprintf('\tInversion Pulse: %d, Set: %d and Inversion Time: %d\n',InversionPulseIndex,MultiBandSetIndex,InversionTimeIndex);
        InversionTimes(InversionTimeIndex,MultiBandSetIndex,InversionPulseIndex)=Header.Time;
    end
    DataIR.Kspace(Kx,Ky,InversionTimeIndex,:,MultiBandSetIndex,InversionPulseIndex)=ComplexData;
end
if (InversionTimeIndex~=DataIR.Info.NumberOfInversionTimes)
    fprintf('jLoadDataIR(): >>>> Warning: Last inversion pulse may not be complete <<<<\n');
end
DataIR.Kspace=permute(DataIR.Kspace,[1 2 3 5 6 4]);

% Find inversion times
TimeStamp=sort(InversionTimes(:));
TimeStamp(find(TimeStamp==0))=[];
if (IsGated)
    DataIR.Info.AverageRR_ms=1e3*mean(diff(TimeStamp));
    HeartRate=60./diff(TimeStamp);
    HeartRate(find(HeartRate==min(HeartRate),1))=[];
    HeartRate(find(HeartRate==max(HeartRate),1))=[];
    DataIR.Info.HeartRateAverage_bpm=mean(HeartRate);
    DataIR.Info.HeartRateStandardDeviation_bpm=std(HeartRate);
    DataIR.Info.MissedBeats=length(find(diff(TimeStamp)>(1.5*mean(diff(TimeStamp)))));
end
DataIR.Info.InversionTimes_ms=1e3*(InversionTimes-repmat(InversionTimes(1,:,:),[DataIR.Info.NumberOfInversionTimes 1 1]))+FirstInversionTime_ms;
DataIR.Info.InversionPulseTimes_ms=squeeze(InversionTimes(1,:,:));
DataIR.Info.InversionPulseTimes_ms=DataIR.Info.InversionPulseTimes_ms-min(DataIR.Info.InversionPulseTimes_ms(:));

% Align data if required
if (exist('AlignKspaceCenter','var')~=1)
    AlignKspaceCenter='CalculateAndApply';
end
if (contains(AlignKspaceCenter,'Calculate','IgnoreCase',true))
    fprintf('jLoadDataIR(): Aligning data\n');
    PhaseRF=DataIR.PhaseRF(:,1,:,:);
    PhaseRF=PhaseRF(:);
    AlignData=reshape(DataIR.Kspace(:,:,1,:,:,:),DataIR.Info.ReadoutSamples,length(PhaseRF),1,DataIR.Info.NumberOfChannels);
    AlignIndex=find(PhaseRF==0);
    N=length(AlignIndex);
    if (N>10)
        ShiftMatrix=repmat(reshape(2i*pi*((1:DataIR.Info.ReadoutSamples)-floor(DataIR.Info.ReadoutSamples/2)-1)/DataIR.Info.ReadoutSamples,DataIR.Info.ReadoutSamples,1),[1 N 1 DataIR.Info.NumberOfChannels]);
        CoilWeight=max(max(abs(AlignData(:,AlignIndex,1,:)),[],1),[],2).^2;
        CoilWeight=CoilWeight./repmat(sum(CoilWeight,4),[1 1 1 DataIR.Info.NumberOfChannels]);
        DataIR.Info.AlignPixel.Offset=linspace(-2,2,100);
        DataIR.Info.AlignPixel.Deviation=zeros(size(DataIR.Info.AlignPixel.Offset),'single');
        for OffsetIndex=1:length(DataIR.Info.AlignPixel.Offset)
            TempData=fftshift(ifft(ifftshift(fftshift(fft(ifftshift(AlignData(:,AlignIndex,1,:),1),[],1),1).*exp(DataIR.Info.AlignPixel.Offset(OffsetIndex)*ShiftMatrix),1),[],1),1);
            PhaseData=exp(1i*angle(TempData(floor(DataIR.Info.ReadoutSamples/2)+1,:,1,:)));
            PhaseDifference=sum(angle(PhaseData.*repmat(conj(sum(PhaseData,2)/N),[1 N 1 1])).^2,2);
            DataIR.Info.AlignPixel.Deviation(OffsetIndex)=sum(PhaseDifference.*CoilWeight);
        end
        DataIR.Info.AlignPixel.AppliedShift=DataIR.Info.AlignPixel.Offset(round(mean(find(DataIR.Info.AlignPixel.Deviation==min(DataIR.Info.AlignPixel.Deviation)))));
        if (contains(AlignKspaceCenter,'Apply','IgnoreCase',true))
            ShiftMatrix=repmat(reshape(2i*pi*((1:DataIR.Info.ReadoutSamples)-floor(DataIR.Info.ReadoutSamples/2)-1)/DataIR.Info.ReadoutSamples,DataIR.Info.ReadoutSamples,1),[1 size(DataIR.Kspace,2) size(DataIR.Kspace,3) size(DataIR.Kspace,4) size(DataIR.Kspace,5) size(DataIR.Kspace,6)]);
            DataIR.Kspace=fftshift(ifft(ifftshift(fftshift(fft(ifftshift(DataIR.Kspace,1),[],1),1).*exp(DataIR.Info.AlignPixel.AppliedShift*ShiftMatrix),1),[],1),1);
        end
    end
end

% Close file if opned in this function
fclose(FileInfo.FileId);


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