% --------------------------------------------------------------------------------
% [KspaceData]=jExportKspace(Data,Theta,Navigator)
% Jason Mendes - UCAIR - University of Utah - 2015
% --------------------------------------------------------------------------------
% Export indices to create kspace.
% --------------------------------------------------------------------------------
function [KspaceData,Kx,Ky]=jExportKspace_v2(Data,Theta,Navigator)

% Find Kspace coordinates
KspaceData=Data;
[Nx,Nr,Nz,Nm,Nc]=size(KspaceData);
Kx=zeros(Nx,Nr,Nz,Nm,'single');
Ky=zeros(Nx,Nr,Nz,Nm,'single');
BaseIndex=repmat(reshape((1:Nx)-floor(Nx/2)-1.0,Nx,1),[1 Nr]);
for MeasIndex=1:Nm
    for IndexZ=1:Nz
        CurrentTheta=repmat(Theta(MeasIndex,:,IndexZ)/180*pi,[Nx 1]);
        Kx(:,:,IndexZ,MeasIndex)=BaseIndex.*cos(CurrentTheta);
        Ky(:,:,IndexZ,MeasIndex)=BaseIndex.*sin(CurrentTheta);
    end
end

% Shift center of Kspace
PhaseMask=repmat(exp(-2i*pi*ndgrid((1:Nx)-floor(Nx/2)-1,1:Nr)/Nx.*0.5),[1 1 Nz Nm Nc]);
KspaceData=fftshift(fft(ifftshift(KspaceData,1),[],1),1).*PhaseMask;
KspaceData=fftshift(ifft(ifftshift(KspaceData,1),[],1),1);

% Fix phase at center of Kspace
for ChannelIndex=1:Nc
    for MeasIndex=1:Nm
        PhaseCorrection=exp(-1i*angle(mean(KspaceData(floor(Nx/2)+1,:,floor(Nz/2)+1,MeasIndex,ChannelIndex),2)));
        KspaceData(:,:,:,MeasIndex,ChannelIndex)=KspaceData(:,:,:,MeasIndex,ChannelIndex).*repmat(PhaseCorrection,[Nx,Nr,Nz,1 1]);
    end
    for IndexZ=1:Nz
        CenterData=KspaceData(floor(Nx/2)+1,:,IndexZ,:,ChannelIndex);
        Mask=1.0*(abs(CenterData)>0);
        MeanValue=sum(CenterData(:).*Mask(:))/sum(Mask(:));
        LocalPhaseCorrection=exp(-1i*angle(CenterData))*exp(1i*angle(MeanValue));    
        KspaceData(:,:,IndexZ,:,ChannelIndex)=KspaceData(:,:,IndexZ,:,ChannelIndex).*repmat(LocalPhaseCorrection,[Nx 1 1 1 1]);
    end
end