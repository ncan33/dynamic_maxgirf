function [imgWindowed, imgCorrected, imgPhase] = phaseSensitiveIR_( images, winSize )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

if nargin<2
    winSize = 4;
end
[imgSize,~,~] = size(images);

xMin = imgSize/2-winSize/2+1;
yMin = imgSize/2-winSize/2+1;

kspace = fftshift2(fft2(fftshift2(images)));
kspaceCtr = kspace(xMin:xMin+winSize-1, yMin:yMin+winSize-1, :);
kspaceLR = padarray(kspaceCtr, [(imgSize-winSize)/2, (imgSize-winSize)/2],0,'both');

imageLR = fftshift2(ifft2(fftshift2(kspaceLR)));

imgPhase = angle(imageLR);
imgCorrected = images ./ exp(1i * imgPhase);
imgWindowed = imgCorrected;
imgWindowed(real(imgCorrected)<0) = 0 + 1i * imag(imgCorrected(real(imgCorrected)<0));

end

