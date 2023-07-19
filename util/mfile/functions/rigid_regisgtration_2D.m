% --------------------------------------------------------------------------------
% [RegisteredData,Offsets]=jRegister2D(UnregisteredData,Mask,ReferenceData)
% Jason Mendes - UCAIR - University of Utah - 2017
% --------------------------------------------------------------------------------
% Use cross correlation to find pixel offsets.
% ----------------------------------------------------------------------
function [RegisteredData]=rigid_regisgtration_2D(UnregisteredData)
Mask = get_mask(abs(UnregisteredData));
% Get data from the selected region
NumberOfIterations=4;
MaxShiftPerIteration=3;
SearchWindow=8;
[Nx,Ny,Nm]=size(UnregisteredData);
KspaceFilter=repmat(reshape(gausswin(Nx,10),Nx,1),[1 Ny]).*repmat(reshape(gausswin(Ny,10),1,Ny),[Nx 1]);
SearchX=(-SearchWindow:SearchWindow)+floor(Nx/2)+1;
Offsets=[];
Offsets.X=zeros(NumberOfIterations,Nm);
Offsets.Y=zeros(NumberOfIterations,Nm);
BaseData = UnregisteredData;
PhaseX=repmat(reshape(-2i*pi*((1:Nx)-floor(Nx/2)-1)/Nx,Nx,1),[1 Ny]);
PhaseY=repmat(reshape(-2i*pi*((1:Ny)-floor(Ny/2)-1)/Ny,1,Ny),[Nx 1]);
RegisteredData=UnregisteredData;
ImageMask=1.0*(conv2(single(Mask),ones(SearchWindow,SearchWindow)/(SearchWindow^2),'same')>0);
ImageMask=conv2(ImageMask,ones(SearchWindow,SearchWindow)/(SearchWindow^2),'same');
for MeasIndex=Nm:-1:1
    for IterIndex=1:NumberOfIterations
        CurrentImage=jFFT(ImageMask.*BaseData(:,:,MeasIndex));
        ReferenceImage=jFFT(ImageMask.*BaseData(:,:,min(MeasIndex+1,Nm)));
        FilteredCrossCorrelation=abs(jIFFT(exp(1i*angle(ReferenceImage.*conj(CurrentImage))).*KspaceFilter));
        ReducedCrossCorrelation=FilteredCrossCorrelation(SearchX,SearchX);
        PeakLocation=find(ReducedCrossCorrelation(:)==max(ReducedCrossCorrelation(:)));
        [IndexX,IndexY]=ind2sub(size(ReducedCrossCorrelation),PeakLocation);
        IndexX=min(max(mean(IndexX)-floor(length(SearchX)/2)+floor(Nx/2),2),Nx-1);
        IndexY=min(max(mean(IndexY)-floor(length(SearchX)/2)+floor(Ny/2),2),Ny-1);
        LocalPeak=sum(FilteredCrossCorrelation((-1:1)+IndexX,(-1:1)+IndexY,:),2);
        Offsets.X(IterIndex,MeasIndex)=max(min(IndexX-0.5*(LocalPeak(3)-LocalPeak(1))/(LocalPeak(1)+LocalPeak(3)-2*LocalPeak(2))-floor(Nx/2)-1,MaxShiftPerIteration),-MaxShiftPerIteration);
        LocalPeak=sum(FilteredCrossCorrelation((-1:1)+IndexX,(-1:1)+IndexY,:),1);
        Offsets.Y(IterIndex,MeasIndex)=max(min(IndexY-0.5*(LocalPeak(3)-LocalPeak(1))/(LocalPeak(1)+LocalPeak(3)-2*LocalPeak(2))-floor(Ny/2)-1,MaxShiftPerIteration),-MaxShiftPerIteration);
        BaseData(:,:,MeasIndex)=abs(jIFFT(jFFT(BaseData(:,:,MeasIndex)).*exp(PhaseX.*Offsets.X(IterIndex,MeasIndex)+PhaseY.*Offsets.Y(IterIndex,MeasIndex))));
        RegisteredData(:,:,MeasIndex)=abs(jIFFT(jFFT(RegisteredData(:,:,MeasIndex)).*exp(PhaseX.*Offsets.X(IterIndex,MeasIndex)+PhaseY.*Offsets.Y(IterIndex,MeasIndex))));   
    end
end

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