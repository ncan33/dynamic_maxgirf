function Image_PD = rigid_reg_PD_3D_yt(Image_PD,Image,Mask)
[sx,sy,sz] = size(Image_PD);
Image = mean(Image,4);

for i=1:sz
    Mask(:,:,i) = bwdist(Mask(:,:,i))<10;
end
Mask = imgaussfilt3(single(Mask),2);
Image_masked = Image.*Mask;
Image_masked_fft = fft3(Image_masked);

PD_fft = fft3(Image_PD.*Mask);

cross_corr = Image_masked_fft.*conj(PD_fft);
cross_corr = ifft3(cross_corr);
cross_corr = fftshift3(cross_corr);
cross_corr = abs(cross_corr);
cross_corr = reshape(cross_corr,[sx*sy*sz,1]);
[~,max_idx] = max(cross_corr,[],1);
[x,y,z] = ind2sub([sx,sy,sz],max_idx);


cross_corr = reshape(cross_corr,[sx,sy,sz]);

LocalPeak = cross_corr(x+(-1:1),y+(-1:1),z+(-1:1));

LocalPeakX = squeeze(sum(sum(LocalPeak,2),3));
LocalPeakY = squeeze(sum(sum(LocalPeak,1),3));
LocalPeakZ = squeeze(sum(sum(LocalPeak,1),2));
OffsetX = x - floor(sx/2) - 1 - 0.5*(LocalPeakX(3)-LocalPeakX(1))./(LocalPeakX(1)+LocalPeakX(3)-2*LocalPeakX(2));
OffsetY = y - floor(sy/2) - 1 - 0.5*(LocalPeakY(3)-LocalPeakY(1))./(LocalPeakY(1)+LocalPeakY(3)-2*LocalPeakY(2));
OffsetZ = z - floor(sz/2) - 1 - 0.5*(LocalPeakZ(3)-LocalPeakZ(1))./(LocalPeakZ(1)+LocalPeakZ(3)-2*LocalPeakZ(2));

PhaseX = reshape(-2i*pi*((1:sx)-floor(sx/2)-1)/sx,[sx,1,1]).*OffsetX;
PhaseY = reshape(-2i*pi*((1:sy)-floor(sy/2)-1)/sy,[1,sy,1]).*OffsetY;
PhaseZ = reshape(-2i*pi*((1:sz+2)-floor((sz+2)/2)-1)/(sz+2),[1,1,sz+2]).*OffsetZ;

Image_PD = cat(3,Image_PD(:,:,1),Image_PD,Image_PD(:,:,end));
Image_PD = fftshift3(ifft3(ifftshift3(fftshift3(fft3(ifftshift3(Image_PD))).*exp(PhaseX+PhaseY+PhaseZ))));
Image_PD = Image_PD(:,:,2:end-1,:);
Image_PD = abs(Image_PD);