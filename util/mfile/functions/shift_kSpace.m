
function kSpace_out = shift_kSpace(kSpace_in,amount)

kSpace_out = kSpace_in;
[sx,nor,nf,Nm,nc]=size(kSpace_out);

% Shift center of Kspace
PhaseMask=repmat(exp(-2i*pi*ndgrid((1:sx)-floor(sx/2)-1,1:nor)/sx.*amount),[1 1 nf Nm nc]);
kSpace_out=fftshift(fft(ifftshift(kSpace_out,1),[],1),1).*PhaseMask;
kSpace_out=fftshift(ifft(ifftshift(kSpace_out,1),[],1),1);

% Fix phase at center of Kspace
for ChannelIndex=1:nc
    for MeasIndex=1:Nm
        PhaseCorrection=exp(-1i*angle(mean(kSpace_out(floor(sx/2)+1,:,floor(nf/2)+1,MeasIndex,ChannelIndex),2)));
        kSpace_out(:,:,:,MeasIndex,ChannelIndex)=kSpace_out(:,:,:,MeasIndex,ChannelIndex).*repmat(PhaseCorrection,[sx,nor,nf,1 1]);
    end
    for IndexZ=1:nf
        CenterData=kSpace_out(floor(sx/2)+1,:,IndexZ,:,ChannelIndex);
        Mask=1.0*(abs(CenterData)>0);
        MeanValue=sum(CenterData(:).*Mask(:))/sum(Mask(:));
        LocalPhaseCorrection=exp(-1i*angle(CenterData))*exp(1i*angle(MeanValue));    
        kSpace_out(:,:,IndexZ,:,ChannelIndex)=kSpace_out(:,:,IndexZ,:,ChannelIndex).*repmat(LocalPhaseCorrection,[sx 1 1 1 1]);
    end
end