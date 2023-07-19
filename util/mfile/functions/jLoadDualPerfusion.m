% --------------------------------------------------------------------------------
% [DataAIF,Data3D,Data2D]=jLoadDualPerfusion(FileName)
% Jason Mendes - UCAIR - University of Utah - 2017
% --------------------------------------------------------------------------------
% Load data from measurement file for interleaved 2D/3D perfusion data.
% --------------------------------------------------------------------------------
function [DataAIF,Data3D,Data2D]=jLoadDualPerfusion(FileName)

% Open file
fprintf('jLoadDualPerfusion(): Opening file %s\n',FileName);
FileInfo=jGetFileInfo(FileName);
DataFormat='Compatible';

% Determine configuration and initialize
SatNormal=0;
SatEfficient=1;
Kernel2D=0;
Kernel3D=1;
KernelAIF=2;
N=length(FileInfo.Protocol);
fprintf('jLoadDualPerfusion(): Loading configuration information\n');
LineCount=zeros(1,3);
TimeStamp=zeros(1,length(FileInfo.Protocol(N).OffsetImage),'single');
MeasurementType=-ones(1,FileInfo.Protocol(N).lRepetitions+1,'single');
MeasurementStartIndex=zeros(1,FileInfo.Protocol(N).lRepetitions+1,'single');
LastRepetitionIndex=-1;
for BlockIndex=1:length(FileInfo.Protocol(N).OffsetImage)
    Header=jGetHeader(FileInfo,BlockIndex);
    TimeStamp(BlockIndex)=Header.Time;
    LineCount(Header.IceParam(8)+1)=LineCount(Header.IceParam(8)+1)+1;
    MeasurementType(Header.Repetition+1)=Header.IceParam(8);
    if (Header.Repetition>LastRepetitionIndex)
        LastRepetitionIndex=Header.Repetition;
        MeasurementStartIndex(Header.Repetition+1)=BlockIndex;
    end
end
TotalMeasurements2D=length(find(MeasurementType==Kernel2D));
TotalMeasurements3D=length(find((MeasurementType==KernelAIF)|(MeasurementType==Kernel3D)));
if ((TotalMeasurements3D+TotalMeasurements2D)~=(FileInfo.Protocol(N).lRepetitions+1))
    fprintf('jLoadDualPerfusion(): Warning: Data file might be corrupt\n')
    fprintf('jLoadDualPerfusion(): Loading in %d of %d measurements\n',TotalMeasurments3D+TotalMeasurments2D,FileInfo.Protocol(N).lRepetitions+1);
end
TriggerTime=TimeStamp(MeasurementStartIndex(1:(TotalMeasurements3D+TotalMeasurements2D)));
TriggerTime=TriggerTime-TriggerTime(1);
try
    TriggerDelay=FileInfo.Protocol(N).alTD(1)/1e3;
catch
    TriggerDelay=0;
end

% Determine 3D setup
Data3D=[];
Data3D.Info.NumberOfSets=FileInfo.Protocol(N).sWipMemBlock.alFree(2);
Data3D.Info.SaturationRecoveryTime=FileInfo.Protocol(N).sWipMemBlock.alFree(3);
Data3D.Info.SetSpacing3D=FileInfo.Protocol(N).sWipMemBlock.alFree(4);
Data3D.Info.ProtonDensityScans=FileInfo.Protocol(N).sWipMemBlock.alFree(9);
Data3D.Info.ProtonDensityFlipAngle=FileInfo.Protocol(N).sWipMemBlock.alFree(10);
Data3D.Info.VariableDensityRatio=FileInfo.Protocol(N).sWipMemBlock.alFree(11)/100;
Data3D.Info.AsymmetricRatio=FileInfo.Protocol(N).sWipMemBlock.alFree(12)/100;
Data3D.Info.RepetitionTime=FileInfo.Protocol(N).alTR(4)/1e3;
Data3D.Info.AcquireTime=FileInfo.Protocol(N).alTR(1)/1e3;
Data3D.Info.EchoTime=FileInfo.Protocol(N).alTE(4)/1e3;
Data3D.Info.LinesPerSet=ceil(LineCount(Kernel3D+1)/TotalMeasurements3D/Data3D.Info.NumberOfSets);
Data3D.Info.NumberOfMeasurements=TotalMeasurements3D;
Data3D.Info.TriggerTime=TriggerTime(find((MeasurementType==KernelAIF)|(MeasurementType==Kernel3D)));
Data3D.Info.AverageRR=mean(diff(TriggerTime));
Data3D.Info.StdDevRR=std(diff(TriggerTime));
Data3D.Info.TriggerDelay=TriggerDelay;
Data3D.Info.Partitions=FileInfo.Protocol(N).sKSpace.lPartitions;
Data3D.Info.LinesInCentralPartition=FileInfo.Protocol(N).sKSpace.lRadialViews;
Data3D.Info.FlipAngle=FileInfo.Protocol(N).adFlipAngleDegree;
Data3D.Info.ReadoutSamples=2*FileInfo.Protocol(N).sKSpace.lBaseResolution;
Data3D.Info.NumberOfChannels=FileInfo.Protocol(N).lNumberOfChannels;
Data3D.Info.FOV=FileInfo.Protocol(N).sSliceArray(1).asSlice(1).dReadoutFOV;
Data3D.Info.SliceThickness=FileInfo.Protocol(N).sSliceArray(1).asSlice(1).dThickness;

% Determine AIF setup
DataAIF.Info.LinesPerSet=FileInfo.Protocol(N).sWipMemBlock.alFree(5);
DataAIF.Info.SliceThickness=FileInfo.Protocol(N).sWipMemBlock.alFree(6);
DataAIF.Info.SaturationRecoveryTime=FileInfo.Protocol(N).sWipMemBlock.alFree(7);
DataAIF.Info.NumberOfContrasts=FileInfo.Protocol(N).sWipMemBlock.alFree(8);
DataAIF.Info.ProtonDensityScans=FileInfo.Protocol(N).sWipMemBlock.alFree(9);
DataAIF.Info.ProtonDensityFlipAngle=FileInfo.Protocol(N).sWipMemBlock.alFree(10);
DataAIF.Info.VariableDensityRatio=FileInfo.Protocol(N).sWipMemBlock.alFree(11)/100;
DataAIF.Info.AsymmetricRatio=FileInfo.Protocol(N).sWipMemBlock.alFree(12)/100;
DataAIF.Info.RepetitionTime=FileInfo.Protocol(N).alTR(2)/1e3;
DataAIF.Info.AcquireTime=FileInfo.Protocol(N).alTR(1)/1e3;
DataAIF.Info.EchoTime=FileInfo.Protocol(N).alTE(2)/1e3;
if (length(FileInfo.Protocol(N).sWipMemBlock.alFree)>12)
    SatType=FileInfo.Protocol(N).sWipMemBlock.alFree(12);
else
    SatType=0;
end
switch SatType
    case SatNormal
        Data.Info.SequenceType='Normal: AIF before 3D and during 3D SRT';
    case SatEfficient
        Data.Info.SequenceType='Compressed: AIF during 3D SRT';
    otherwise
        Data.Info.SequenceType='Unknown';
end
DataAIF.Info.NumberOfSets=ceil(LineCount(KernelAIF+1)/TotalMeasurements3D/DataAIF.Info.LinesPerSet/DataAIF.Info.NumberOfContrasts);
DataAIF.Info.NumberOfMeasurements=TotalMeasurements3D;
DataAIF.Info.TriggerTime=TriggerTime(find((MeasurementType==KernelAIF)|(MeasurementType==Kernel3D)));
DataAIF.Info.AverageRR=mean(diff(TriggerTime));
DataAIF.Info.StdDevRR=std(diff(TriggerTime));
DataAIF.Info.TriggerDelay=TriggerDelay;
DataAIF.Info.FlipAngle=FileInfo.Protocol(N).adFlipAngleDegree;
DataAIF.Info.ReadoutSamples=2*FileInfo.Protocol(N).sKSpace.lBaseResolution;
DataAIF.Info.NumberOfChannels=FileInfo.Protocol(N).lNumberOfChannels;
DataAIF.Info.FOV=FileInfo.Protocol(N).sSliceArray(1).asSlice(1).dReadoutFOV;
if (SatType==SatNormal)
    DataAIF.Info.FileOffset=[0 DataAIF.Info.LinesPerSet*DataAIF.Info.NumberOfContrasts 2*DataAIF.Info.LinesPerSet*DataAIF.Info.NumberOfContrasts+Data3D.Info.LinesPerSet];
    Data3D.Info.FileOffset=[2*DataAIF.Info.LinesPerSet*DataAIF.Info.NumberOfContrasts 3*DataAIF.Info.LinesPerSet*DataAIF.Info.NumberOfContrasts+Data3D.Info.LinesPerSet];
else
    DataAIF.Info.FileOffset=[0 DataAIF.Info.LinesPerSet*DataAIF.Info.NumberOfContrasts+Data3D.Info.LinesPerSet];
    Data3D.Info.FileOffset=[DataAIF.Info.LinesPerSet*DataAIF.Info.NumberOfContrasts 2*DataAIF.Info.LinesPerSet*DataAIF.Info.NumberOfContrasts+Data3D.Info.LinesPerSet];
end

% Determine 2D setup
Data2D.Info.LinesPerSet=ceil(LineCount(Kernel2D+1)/TotalMeasurements2D);
Data2D.Info.SliceThickness=FileInfo.Protocol(N).sWipMemBlock.alFree(6);
Data2D.Info.SaturationRecoveryTime=FileInfo.Protocol(N).sWipMemBlock.alFree(7);
Data2D.Info.ProtonDensityScans=FileInfo.Protocol(N).sWipMemBlock.alFree(9);
Data2D.Info.ProtonDensityFlipAngle=FileInfo.Protocol(N).sWipMemBlock.alFree(10);
Data2D.Info.VariableDensityRatio=FileInfo.Protocol(N).sWipMemBlock.alFree(11)/100;
Data2D.Info.AsymmetricRatio=FileInfo.Protocol(N).sWipMemBlock.alFree(12)/100;
Data2D.Info.RepetitionTime=FileInfo.Protocol(N).alTR(3)/1e3;
Data2D.Info.AcquireTime=FileInfo.Protocol(N).alTR(1)/1e3;
Data2D.Info.EchoTime=FileInfo.Protocol(N).alTE(3)/1e3;
Data2D.Info.NumberOfMeasurements=TotalMeasurements3D;
Data2D.Info.TriggerTime=TriggerTime(find((MeasurementType==KernelAIF)|(MeasurementType==Kernel3D)));
Data2D.Info.AverageRR=mean(diff(TriggerTime));
Data2D.Info.StdDevRR=std(diff(TriggerTime));
Data2D.Info.TriggerDelay=TriggerDelay;
Data2D.Info.FlipAngle=FileInfo.Protocol(N).adFlipAngleDegree;
Data2D.Info.ReadoutSamples=2*FileInfo.Protocol(N).sKSpace.lBaseResolution;
Data2D.Info.NumberOfChannels=FileInfo.Protocol(N).lNumberOfChannels;
Data2D.Info.FOV=FileInfo.Protocol(N).sSliceArray(1).asSlice(1).dReadoutFOV;
Data2D.Info.Multiband=3;
Data2D.Info.NumberOfSlices=3;
Data2D.Info.SliceDistance=0.5*FileInfo.Protocol(N).sSliceArray(1).asSlice(1).dThickness;

% Read in all AIF data
CurrentMeasurements=find((MeasurementType==KernelAIF)|(MeasurementType==Kernel3D));
for SetIndex=1:DataAIF.Info.NumberOfSets
    if (strcmpi(DataFormat,'Compatible'))
        Kspace=complex(0,zeros(DataAIF.Info.ReadoutSamples,DataAIF.Info.LinesPerSet,1,DataAIF.Info.NumberOfMeasurements,DataAIF.Info.NumberOfContrasts,DataAIF.Info.NumberOfChannels,'single'));
        Theta=zeros(DataAIF.Info.NumberOfMeasurements,DataAIF.Info.LinesPerSet,'single');
    else
        Kspace=complex(0,zeros(DataAIF.Info.ReadoutSamples,DataAIF.Info.LinesPerSet,DataAIF.Info.NumberOfMeasurements,DataAIF.Info.NumberOfChannels,DataAIF.Info.NumberOfContrasts,'single'));
        Theta=zeros(DataAIF.Info.LinesPerSet,1,DataAIF.Info.NumberOfMeasurements,'single');
    end
    for MeasIndex=1:DataAIF.Info.NumberOfMeasurements
        fprintf('Reading in AIF from measurement %d of %d (Set %d of %d)\n',MeasIndex,DataAIF.Info.NumberOfMeasurements,SetIndex,DataAIF.Info.NumberOfSets);
        DataOffset=MeasurementStartIndex(CurrentMeasurements(MeasIndex))+DataAIF.Info.FileOffset(SetIndex);
        for LineIndex=1:(DataAIF.Info.LinesPerSet*DataAIF.Info.NumberOfContrasts)
            [ComplexData,Header]=jGetDataAndHeader(FileInfo,DataOffset+LineIndex-1);
            if (Header.IceParam(8)~=KernelAIF)
                fprintf('jLoadDualPerfusion(): AIF line tagged wrong\n');
            end
            Ky=Header.Line+1;
            Kx=(1:Header.Samples)-Header.Samples+DataAIF.Info.ReadoutSamples;
            Header.Echo=rem(LineIndex-1,DataAIF.Info.NumberOfContrasts);
            if (strcmpi(DataFormat,'Compatible'))
                Theta(MeasIndex,Ky)=Header.IceParam(5)/100;
                if (rem(Header.Echo,2))
                    Kspace(Kx,Ky,1,MeasIndex,Header.Echo+1,:)=reshape(ComplexData(Header.Samples:-1:1,:,:,:),Header.Samples,1,1,1,1,DataAIF.Info.NumberOfChannels);
                else
                    Kspace(Kx,Ky,1,MeasIndex,Header.Echo+1,:)=reshape(ComplexData,Header.Samples,1,1,1,1,DataAIF.Info.NumberOfChannels);
                end
            else
                Theta(Ky,1,MeasIndex)=Header.IceParam(5)/100;
                if (rem(Header.Echo,2))
                    Kspace(Kx,Ky,MeasIndex,:,Header.Echo+1)=ComplexData(Header.Samples:-1:1,:,:,:);
                else
                    Kspace(Kx,Ky,MeasIndex,:,Header.Echo+1)=ComplexData;
                end
            end
        end
    end
    switch (SetIndex)
        case 1
            if (SatType==SatNormal)
                if (strcmpi(DataFormat,'Compatible')) 
                    DataAIF.AIF_Normal=Kspace;
                    DataAIF.AIF_Normal_Theta=Theta;
                else
                    DataAIF.Normal.Kspace=Kspace;
                    DataAIF.Normal.Theta=Theta;
                end
            else
                if (strcmpi(DataFormat,'Compatible')) 
                    DataAIF.AIF_Systole=Kspace;
                    DataAIF.AIF_Systole_Theta=Theta;
                else
                    DataAIF.Systole.Kspace=Kspace;
                    DataAIF.Systole.Theta=Theta;
                end
            end
        case 2
            if (SatType==SatNormal)
                if (strcmpi(DataFormat,'Compatible')) 
                    DataAIF.AIF_Systole=Kspace;
                    DataAIF.AIF_Systole_Theta=Theta;
                else
                    DataAIF.Systole.Kspace=Kspace;
                    DataAIF.Systole.Theta=Theta;
                end
            else
                if (strcmpi(DataFormat,'Compatible')) 
                    DataAIF.AIF_Diastole=Kspace;
                    DataAIF.AIF_Diastole_Theta=Theta;
                else
                    DataAIF.Diastole.Kspace=Kspace;
                    DataAIF.Diastole.Theta=Theta;
                end
            end
        otherwise
            if (strcmpi(DataFormat,'Compatible')) 
                DataAIF.AIF_Diastole=Kspace;
                DataAIF.AIF_Diastole_Theta=Theta;
            else
                DataAIF.Diastole.Kspace=Kspace;
                DataAIF.Diastole.Theta=Theta;
            end
    end
end

% Read in all 3D data
CurrentMeasurements=find((MeasurementType==KernelAIF)|(MeasurementType==Kernel3D));
for SetIndex=1:Data3D.Info.NumberOfSets
    if (strcmpi(DataFormat,'Compatible'))
        Kspace=complex(0,zeros(Data3D.Info.ReadoutSamples,Data3D.Info.LinesInCentralPartition,Data3D.Info.Partitions,Data3D.Info.NumberOfMeasurements,Data3D.Info.NumberOfChannels,'single'));
        Theta=zeros(Data3D.Info.NumberOfMeasurements,Data3D.Info.LinesInCentralPartition,Data3D.Info.Partitions,'single');
    else
        Kspace=complex(0,zeros(Data3D.Info.ReadoutSamples,Data3D.Info.LinesInSet,Data3D.Info.NumberOfMeasurements,Data3D.Info.NumberOfChannels,'single'));
        Theta=zeros(Data3D.Info.LinesPerSet,Data3D.Info.NumberOfMeasurements,'single');
        Kz=zeros(Data3D.Info.LinesPerSet,Data3D.Info.NumberOfMeasurements,'single');
    end
    for MeasIndex=1:Data3D.Info.NumberOfMeasurements
        fprintf('Reading in 3D from measurement %d of %d (Set %d of %d)\n',MeasIndex,Data3D.Info.NumberOfMeasurements,SetIndex,Data3D.Info.NumberOfSets);
        DataOffset=MeasurementStartIndex(CurrentMeasurements(MeasIndex))+Data3D.Info.FileOffset(SetIndex);
        for LineIndex=1:Data3D.Info.LinesPerSet
            [ComplexData,Header]=jGetDataAndHeader(FileInfo,DataOffset+LineIndex-1);
            if (Header.IceParam(8)~=Kernel3D)
                fprintf('jLoadDualPerfusion(): 3D line tagged wrong\n');
            end
            Ky=Header.Line+1;
            Kx=(1:Header.Samples)-Header.Samples+DataAIF.Info.ReadoutSamples;
            if (strcmpi(DataFormat,'Compatible'))
                Kz=Header.Partition-Header.KspacePartitionCenter+floor(Data3D.Info.Partitions/2)+1;
                Theta(MeasIndex,Ky,Kz)=Header.IceParam(5)/100;
                Kspace(Kx,Ky,Kz,MeasIndex,:)=reshape(ComplexData,Header.Samples,1,1,1,Data3D.Info.NumberOfChannels);
            else
                Kz(LineIndex,MeasIndex)=Header.Partition-Header.KspacePartitionCenter+floor(Data3D.Info.Partitions/2)+1;
                Theta(LineIndex,MeasIndex)=Header.IceParam(5)/100;
                Kspace(Kx,LineIndex,MeasIndex,:)=ComplexData;
            end
        end
    end
    switch (SetIndex)
        case 1
            if (strcmpi(DataFormat,'Compatible')) 
                Data3D.Systole=Kspace;
                Data3D.Systole_Theta=Theta;
            else
                Data3D.Systole.Kspace=Kspace;
                Data3D.Systole.Theta=Theta;
                Data3D.Systole.Kz=Theta;
            end
        otherwise
            if (strcmpi(DataFormat,'Compatible')) 
                Data3D.Diastole=Kspace;
                Data3D.Diastole_Theta=Theta;
            else
                Data3D.Diastole.Kspace=Kspace;
                Data3D.Diastole.Theta=Theta;
                Data3D.Diastole.Kz=Theta;
            end
    end
end

% Read in all 2D data
CurrentMeasurements=find((MeasurementType==Kernel2D));
Data2D.Kspace=complex(0,zeros(Data2D.Info.ReadoutSamples,Data2D.Info.LinesPerSet,Data2D.Info.NumberOfMeasurements,Data2D.Info.NumberOfChannels,'single'));
Data2D.Theta=zeros(Data2D.Info.LinesPerSet,Data2D.Info.NumberOfMeasurements,'single');
Data2D.PhaseModulation=zeros(Data2D.Info.LinesPerSet,Data2D.Info.NumberOfMeasurements,'single');
for MeasIndex=1:Data2D.Info.NumberOfMeasurements
    fprintf('Reading in 2D from measurement %d of %d\n',MeasIndex,Data2D.Info.NumberOfMeasurements);
    DataOffset=MeasurementStartIndex(CurrentMeasurements(MeasIndex));
    for LineIndex=1:Data2D.Info.LinesPerSet
        [ComplexData,Header]=jGetDataAndHeader(FileInfo,DataOffset+LineIndex-1);
        if (Header.IceParam(8)~=Kernel2D)
            fprintf('jLoadDualPerfusion(): 2D line tagged wrong\n');
        end
        Kx=(1:Header.Samples)-Header.Samples+DataAIF.Info.ReadoutSamples;
        Data2D.Theta(LineIndex,MeasIndex)=Header.IceParam(5)/100;
        Data2D.PhaseModulation(LineIndex,MeasIndex)=Header.IceParam(7);
        Data2D.Kspace(Kx,LineIndex,MeasIndex,:)=ComplexData;
    end
end
if (strcmpi(DataFormat,'Compatible'))
    Data2D.Theta=permute(Data2D.Theta,[2 1]);
    Data2D.PhaseModulation=permute(Data2D.PhaseModulation,[2 1]);
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