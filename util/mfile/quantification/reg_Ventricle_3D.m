function Image_reg = reg_Ventricle_3D(Image,Mask_all)

[sx,sy,sz,nof] = size(Image);
curve = Image(repmat(Mask_all,[1,1,1,nof]));
Npoints = length(curve)/nof;
curve = reshape(curve,[Npoints,nof]);
max_int = max(curve,[],2);
curve_mean = mean(curve,1);

[~,peak_enhan_frame] = max(smooth(curve_mean));
[~,enhan_begin] = max(abs(diff(smooth(curve_mean))));
N_pre = round(enhan_begin*0.8);

curve_mean = curve_mean - mean(curve_mean(1:N_pre));
% AIF_param = Quant.fit_AIF(1:nof,double(curve_mean));
% AIF_model = Quant.AIF_model(AIF_param,1:nof);

model_image = zeros(size(Image));
pre_int = mean(curve(:,1:N_pre),2);
%AIF_model_all = AIF_model/max(AIF_model).*(max_int-pre_int) + pre_int;
curve_all = curve_mean/max(curve_mean).*(max_int-pre_int) + pre_int;
model_image(repmat(Mask_all,[1,1,1,nof])) = curve_all(:);

for i=1:sz
    Mask(:,:,i) = bwdist(Mask_all(:,:,i))<=5;
end
Mask = imgaussfilt3(single(Mask),2);
Image_masked = Image.*Mask;
Image_masked_fft = fft3(Image_masked);

Model_fft = fft3(model_image);


cross_corr = Image_masked_fft.*conj(Model_fft);
cross_corr = ifft3(cross_corr);
cross_corr = fftshift3(cross_corr);
cross_corr = abs(cross_corr);
cross_corr = reshape(cross_corr,[sx*sy*sz,nof]);
[~,max_idx] = max(cross_corr,[],1);
[x,y,z] = ind2sub([sx,sy,sz],max_idx);


cross_corr = reshape(cross_corr,[sx,sy,sz,nof]);
for i=1:nof
    LocalPeak(:,:,:,i) = cross_corr(x(i)+(-1:1),y(i)+(-1:1),z(i)+(-1:1),i);
end
LocalPeakX = squeeze(sum(sum(LocalPeak,2),3));
LocalPeakY = squeeze(sum(sum(LocalPeak,1),3));
LocalPeakZ = squeeze(sum(sum(LocalPeak,1),2));
OffsetX = x - floor(sx/2) - 1 - 0.5*(LocalPeakX(3,:)-LocalPeakX(1,:))./(LocalPeakX(1,:)+LocalPeakX(3,:)-2*LocalPeakX(2,:));
OffsetY = y - floor(sy/2) - 1 - 0.5*(LocalPeakY(3,:)-LocalPeakY(1,:))./(LocalPeakY(1,:)+LocalPeakY(3,:)-2*LocalPeakY(2,:));
OffsetZ = z - floor(sz/2) - 1 - 0.5*(LocalPeakZ(3,:)-LocalPeakZ(1,:))./(LocalPeakZ(1,:)+LocalPeakZ(3,:)-2*LocalPeakZ(2,:));

PhaseX = reshape(-2i*pi*((1:sx)-floor(sx/2)-1)/sx,[sx,1,1]).*reshape(OffsetX,[1,1,1,nof]);
PhaseY = reshape(-2i*pi*((1:sy)-floor(sy/2)-1)/sy,[1,sy,1]).*reshape(OffsetY,[1,1,1,nof]);
PhaseZ = reshape(-2i*pi*((1:sz+2)-floor((sz+2)/2)-1)/(sz+2),[1,1,sz+2]).*reshape(OffsetZ,[1,1,1,nof]);

Image = cat(3,Image(:,:,1,:),Image,Image(:,:,end,:));
Image_reg = fftshift3(ifft3(ifftshift3(fftshift3(fft3(ifftshift3(Image))).*exp(-(PhaseX+PhaseY+PhaseZ)))));
Image_reg = Image_reg(:,:,2:end-1,:);
Image_reg = abs(Image_reg);