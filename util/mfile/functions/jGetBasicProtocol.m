% --------------------------------------------------------------------------------
%  ProtocolData=jGetBasicProtocol(BasicData,Sequence)
% Jason Mendes - UCAIR
% --------------------------------------------------------------------------------
% Define the basic protocol
% --------------------------------------------------------------------------------
function ProtocolData=jGetBasicProtocol(BasicData,Sequence)

% Get all common protocol information
ProtocolData=[];
ProtocolData.Sequence=Sequence;
ProtocolData.SiemensHeader=BasicData.Protocol(length(BasicData.Protocol));
ProtocolData.TotalNumberOfLines=length(ProtocolData.SiemensHeader.OffsetImage);
ProtocolData.NumberOfMeasurements=ProtocolData.SiemensHeader.lRepetitions+1;
if (ProtocolData.SiemensHeader.sTXSPEC.ucRFPulseType==1)
    ProtocolData.PulseType='Fast RF';
elseif (ProtocolData.SiemensHeader.sTXSPEC.ucRFPulseType==4)
    ProtocolData.PulseType='Low SAR RF';
else
    ProtocolData.PulseType='Normal RF';
end
ProtocolData.ReadoutSamples=2*ProtocolData.SiemensHeader.sKSpace.lBaseResolution;
ProtocolData.NumberOfChannels=ProtocolData.SiemensHeader.lNumberOfChannels;
ProtocolData.FlipAngle=ProtocolData.SiemensHeader.adFlipAngleDegree;
ProtocolData.FOV_mm=ProtocolData.SiemensHeader.sSliceArray(1).asSlice(1).dReadoutFOV;
ProtocolData.RadialViews=ProtocolData.SiemensHeader.sKSpace.lRadialViews;
try 
    ProtocolData.ProtonDensityScans=ProtocolData.SiemensHeader.lProtonDensMap;
catch
    ProtocolData.ProtonDensityScans=0;
end
if (ProtocolData.ProtonDensityScans>0)
    ProtocolData.ProtonDensityFlipAngle=2;
end
ParameterWIP=zeros(1,14);
ParameterWIP(1:(length(ProtocolData.SiemensHeader.sWipMemBlock.alFree)-1))=ProtocolData.SiemensHeader.sWipMemBlock.alFree(2:length(ProtocolData.SiemensHeader.sWipMemBlock.alFree));
ProtocolData.SiemensHeader.sWipMemBlock.alFree=ParameterWIP;
TR=zeros(1,20);
TR(1:min(length(ProtocolData.SiemensHeader.alTR),3))=ProtocolData.SiemensHeader.alTR(1:min(length(ProtocolData.SiemensHeader.alTR),3))*1e-3;
TE=zeros(1,20);
TE(1:min(length(ProtocolData.SiemensHeader.alTE),3))=ProtocolData.SiemensHeader.alTE(1:min(length(ProtocolData.SiemensHeader.alTE),3))*1e-3;
ProtocolData.TimePerMeasurement_ms=TR(1);
ProtocolData.SiemensHeader.alTR=TR;
ProtocolData.SiemensHeader.alTE=TE;

% Get physio information
if (ProtocolData.SiemensHeader.sPhysioImaging.lSignal1==2)
    FirstLineOfMeasurement=[0 find(diff(BasicData.Measurement)>0)]+1;
    TimeStamps=diff(BasicData.TimeStamp(FirstLineOfMeasurement)-BasicData.PMU(FirstLineOfMeasurement));
    HeartRate=60./TimeStamps;
    HeartRate(find(HeartRate==min(HeartRate),1))=[];
    HeartRate(find(HeartRate==max(HeartRate),1))=[];
    ProtocolData.HeartRate_bpm=round(mean(HeartRate));
    ProtocolData.HeartRate_Deviation_bpm=round(10*std(HeartRate))/10;
    ProtocolData.MissedTriggers=sum(HeartRate>(1.5*ProtocolData.HeartRate_bpm));
    try
        ProtocolData.TriggerDelay_ms=Data.Protocol.alTD(1)*1e-3;
    catch
        ProtocolData.TriggerDelay_ms=[];
    end    
end
