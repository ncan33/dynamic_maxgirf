% --------------------------------------------------------------------------------
% [KspaceData, Kx, Ky] = yExportKspace(Data,Theta,Navigator)
% Ye Tian - UCAIR - University of Utah - 2017
% --------------------------------------------------------------------------------
% Export indices to create kspace.
% --------------------------------------------------------------------------------
function [KspaceData,Kx,Ky] = yExportKspace(Data,Theta,center_shift)

% Find Kspace coordinates
KspaceData = Data;
[Nx,Nr,Nz,Nm,Nc] = size(KspaceData);
Theta = permute(Theta,[4,2,3,1]);
R = (1:Nx)-floor(Nx/2)-1;
R = R.';
Kx = R.*cos(Theta/180*pi);
Ky = R.*sin(Theta/180*pi);

% Shift center of Kspace
Mask = logical(abs(KspaceData));
PhaseMask = exp(-2i*pi*R/Nx*center_shift);
KspaceData = fftshift(fft(ifftshift(KspaceData,1),[],1),1).*PhaseMask;
KspaceData = fftshift(ifft(ifftshift(KspaceData,1),[],1),1);
KspaceData = KspaceData.*Mask;

% Fix phase at center of Kspace
% for ChannelIndex=1:Nc
%     for MeasIndex=1:Nm
%         PhaseCorrection=exp(-1i*angle(mean(KspaceData(floor(Nx/2)+1,:,floor(Nz/2)+1,MeasIndex,ChannelIndex),2)));
%         KspaceData(:,:,:,MeasIndex,ChannelIndex)=KspaceData(:,:,:,MeasIndex,ChannelIndex).*repmat(PhaseCorrection,[Nx,Nr,Nz,1 1]);
%     end
%     for IndexZ=1:Nz
%         CenterData=KspaceData(floor(Nx/2)+1,:,IndexZ,:,ChannelIndex);
%         Mask=1.0*(abs(CenterData)>0);
%         MeanValue=sum(CenterData(:).*Mask(:))/sum(Mask(:));
%         LocalPhaseCorrection=exp(-1i*angle(CenterData))*exp(1i*angle(MeanValue));    
%         KspaceData(:,:,IndexZ,:,ChannelIndex)=KspaceData(:,:,IndexZ,:,ChannelIndex).*repmat(LocalPhaseCorrection,[Nx 1 1 1 1]);
%     end
% end