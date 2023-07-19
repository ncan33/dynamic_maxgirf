function [kSpace,BaseData,Offsets] = rigid_regisgtration_3D(UnregisteredImage,Data)

% Get data from the selected region
Mask = get_mask(sum(abs(UnregisteredImage(:,:,5,:)),4));
if nargin<2
    flag = 0;
else
    flag = 1;
end
if flag
    kSpace = fftshift3(Data.kSpace);
    kSpace = ifft3(kSpace);
    kSpace = fftshift3(kSpace);
    kSpace = fft3(kSpace);
end
NumberOfIterations=4;
MaxShiftPerIteration=3;
SearchWindowX=8;
SearchWindowZ=1;
[Nx,Ny,Nz,Nm]=size(UnregisteredImage);
KspaceFilter=repmat(reshape(gausswin(Nx,10),Nx,1),[1 Ny Nz]).*repmat(reshape(gausswin(Ny,10),1,Ny),[Nx 1 Nz]).*repmat(reshape(gausswin(Nz,10),1,1,Nz),[Nx Ny 1]);
SearchX=(-SearchWindowX:SearchWindowX)+floor(Nx/2)+1;
SearchZ=(-SearchWindowZ:SearchWindowZ)+floor(Nz/2)+1;
Offsets=[];
Offsets.X=zeros(NumberOfIterations,Nm);
Offsets.Y=zeros(NumberOfIterations,Nm);
Offsets.Z=zeros(NumberOfIterations,Nm);
BaseData=UnregisteredImage;
PhaseX=repmat(reshape(-2i*pi*((1:Nx)-floor(Nx/2)-1)/Nx,Nx,1),[1 Ny Nz]);
PhaseY=repmat(reshape(-2i*pi*((1:Ny)-floor(Ny/2)-1)/Ny,1,Ny),[Nx 1 Nz]);
PhaseZ=repmat(reshape(-2i*pi*((1:Nz)-floor(Nz/2)-1)/Nz,1,1,Nz),[Nx Ny 1]);
%RegisteredData=UnregisteredImage;
ImageMaskX=1.0*(conv2(single(Mask),ones(SearchWindowX,SearchWindowX)/(SearchWindowX^2),'same')>0);
ImageMaskX=conv2(ImageMaskX,ones(SearchWindowX,SearchWindowX)/(SearchWindowX^2),'same');
ImageMaskZ=triang(Nz);
ImageMaskZ=min(min(ImageMaskZ(SearchZ)),ImageMaskZ);
ImageMaskZ=reshape(ImageMaskZ/max(ImageMaskZ),1,1,Nz);
ImageMask=repmat(ImageMaskZ,[Nx Ny 1]).*repmat(ImageMaskX,[1 1 Nz]);
for MeasIndex=Nm:-1:1
    for IterIndex=1:NumberOfIterations
        CurrentImage=jFFT(ImageMask.*BaseData(:,:,:,MeasIndex));
        ReferenceImage=jFFT(ImageMask.*BaseData(:,:,:,min(MeasIndex+1,Nm)));
        FilteredCrossCorrelation=abs(jIFFT(exp(1i*angle(ReferenceImage.*conj(CurrentImage))).*KspaceFilter));
        ReducedCrossCorrelation=FilteredCrossCorrelation(SearchX,SearchX,SearchZ);
        PeakLocation=find(ReducedCrossCorrelation(:)==max(ReducedCrossCorrelation(:)),1);
        [IndexX,IndexY,IndexZ]=ind2sub(size(ReducedCrossCorrelation),PeakLocation);
        IndexX=min(max(mean(IndexX)-floor(length(SearchX)/2)+floor(Nx/2),2),Nx-1);
        IndexY=min(max(mean(IndexY)-floor(length(SearchX)/2)+floor(Ny/2),2),Ny-1);
        IndexZ=min(max(mean(IndexZ)-floor(length(SearchZ)/2)+floor(Nz/2),2),Nz-1);
        LocalPeak=sum(sum(FilteredCrossCorrelation((-1:1)+IndexX,(-1:1)+IndexY,(-1:1)+IndexZ),2),3);
        Offsets.X(IterIndex,MeasIndex)=max(min(IndexX-0.5*(LocalPeak(3)-LocalPeak(1))/(LocalPeak(1)+LocalPeak(3)-2*LocalPeak(2))-floor(Nx/2)-1,MaxShiftPerIteration),-MaxShiftPerIteration);
        LocalPeak=sum(sum(FilteredCrossCorrelation((-1:1)+IndexX,(-1:1)+IndexY,(-1:1)+IndexZ),1),3);
        Offsets.Y(IterIndex,MeasIndex)=max(min(IndexY-0.5*(LocalPeak(3)-LocalPeak(1))/(LocalPeak(1)+LocalPeak(3)-2*LocalPeak(2))-floor(Ny/2)-1,MaxShiftPerIteration),-MaxShiftPerIteration);
        LocalPeak=sum(sum(FilteredCrossCorrelation((-1:1)+IndexX,(-1:1)+IndexY,(-1:1)+IndexZ),1),2);
        Offsets.Z(IterIndex,MeasIndex)=max(min(IndexZ-0.5*(LocalPeak(3)-LocalPeak(1))/(LocalPeak(1)+LocalPeak(3)-2*LocalPeak(2))-floor(Nz/2)-1,MaxShiftPerIteration),-MaxShiftPerIteration);       
        BaseData(:,:,:,MeasIndex)=abs(jIFFT(jFFT(BaseData(:,:,:,MeasIndex)).*exp(PhaseX.*Offsets.X(IterIndex,MeasIndex)+PhaseY.*Offsets.Y(IterIndex,MeasIndex)+PhaseZ.*Offsets.Z(IterIndex,MeasIndex))));
        if flag
            kSpace(:,:,:,MeasIndex) = kSpace(:,:,:,MeasIndex).*exp(PhaseX.*Offsets.X(IterIndex,MeasIndex)+PhaseY.*Offsets.Y(IterIndex,MeasIndex)+PhaseZ.*Offsets.Z(IterIndex,MeasIndex));
        end
        %RegisteredData(:,:,:,MeasIndex)=abs(jIFFT(jFFT(RegisteredData(:,:,:,MeasIndex)).*exp(PhaseX.*Offsets.X(IterIndex,MeasIndex)+PhaseY.*Offsets.Y(IterIndex,MeasIndex)+PhaseZ.*Offsets.Z(IterIndex,MeasIndex))));
    end
end
if flag == 0
    kSpace = 0;
    return
end
kSpace = fftshift3(kSpace);
kSpace = ifft3(kSpace);
kSpace = fftshift3(kSpace);
kSpace = fft3(kSpace);
kSpace = kSpace.*Data.mask;

% --------------------------------------------------------------------------------
% F=jFFT(f)
% Jason Mendes - UCAIR - University of Utah - 2012
% --------------------------------------------------------------------------------
% Symmetric FFT.
% --------------------------------------------------------------------------------
function F=jFFT(f)

if (ndims(f)>3)
   error('jFFT(): Data must be 1D, 2D or 3D');
elseif (ndims(f)==3) % 3D
   F=fftshift(fft(fft(fft(ifftshift(f),[],1),[],2),[],3))/sqrt(numel(f));
elseif (min(size(f))>1) % 2D
   F=fftshift(fft2(ifftshift(f)))/sqrt(numel(f));
else % 1D
   F=fftshift(fft(ifftshift(f)))/sqrt(numel(f));
end


% --------------------------------------------------------------------------------
% f=jIFFT(F)
% Jason Mendes - UCAIR - University of Utah - 2012
% --------------------------------------------------------------------------------
% Symmetric inverse FFT.
% --------------------------------------------------------------------------------
function f=jIFFT(F)

if (ndims(F)>3)
   error('jFFT(): Data must be 1D, 2D or 3D');
elseif (ndims(F)==3) % 3D
   f=fftshift(ifft(ifft(ifft(ifftshift(F),[],1),[],2),[],3))*sqrt(numel(F));
elseif (min(size(F))>1) % 2D
   f=fftshift(ifft2(ifftshift(F)))*sqrt(numel(F));
else % 1D
   f=fftshift(ifft(ifftshift(F)))*sqrt(numel(F));
end