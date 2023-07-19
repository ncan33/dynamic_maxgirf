function Mask_MYO = seg_Myo_3D_MB(Image,Mask_LV,threshold)



[sx,sy,sz,nof] = size(Image);
curves = Image(repmat(Mask_LV,[1,1,1,nof]));
Npoints = length(curves)/nof;

curves = reshape(curves,[Npoints,nof]);
curve_mean = mean(curves);
[~,peak_enhan_frame] = max(smooth(curve_mean));
[~,enhan_begin] = max(abs(diff(smooth(curve_mean))));
N_pre = round(enhan_begin*0.8);

curve_mean = curve_mean-mean(curve_mean(1:N_pre));
AIF_param = Quant.fit_AIF(1:nof,double(curve_mean));
AIF_model = Quant.AIF_model(AIF_param,1:nof)';

%% initialize mask myo
Mask_MYO = false(sx,sy,sz);
for i=1:sz
    Mask_MYO(:,:,i) = bwdist(Mask_LV(:,:,i))<=10 & bwdist(Mask_LV(:,:,i))>2;
end
Mask_MYO = Mask_MYO&~Mask_LV;

curves = Image(repmat(Mask_MYO,[1,1,1,nof]));
Npoints = length(curves)/nof;
curves = reshape(curves,[Npoints,nof]);

Nbins = 10;
bins = kmeans(curves,Nbins);
clear curve_mean
for i=1:Nbins
    curve_mean(i,:) = mean(curves(bins==i,:));
end

curve_mean = curve_mean - mean(curve_mean(:,1:N_pre),2);
idx_keep = true(1,Nbins);
idx_keep(min(curve_mean(:,N_pre:end),[],2)<0) = false;
idx_keep(max(curve_mean,[],2)==curve_mean(:,peak_enhan_frame)) = false;
idx_keep(std(curve_mean(:,1:N_pre),1,2)>5) = false;
idx_drop = ~idx_keep;
idx = find(Mask_MYO);
for i=find(idx_drop)
    Mask_MYO(idx(bins==i)) = false;
end
curves = Image(repmat(Mask_MYO,[1,1,1,nof]));
Npoints = length(curves)/nof;
curves = reshape(curves,[Npoints,nof]);
curve_mean = mean(curves,1);
curve_mean = curve_mean - mean(curve_mean(1:N_pre));
pinv_curve = pinv(curve_mean)';
%% refine mask myo

Niter = 10;
for iter=1:Niter
    Mask_analysis = bwdist(Mask_MYO)<=3 & ~Mask_MYO;
    idx = find(Mask_analysis);
    [x,y,z] = ind2sub([sx,sy,sz],idx);
    N_analysis = length(idx);
    pd = zeros(1,N_analysis);
    for i=1:N_analysis
        curve_temp = squeeze(Image(x(i),y(i),z(i),:));
        curve_temp = curve_temp - mean(curve_temp(1:N_pre));
        scale = 1/(pinv_curve*curve_temp);
        pd(i) = mean(abs(curve_temp'*scale - curve_mean)./max(curve_mean));
    end
    keep = pd<threshold;
    Mask_keep = false(sx,sy,sz);
    Mask_keep(idx(keep)) = true;
    Mask_MYO = (Mask_MYO & ~Mask_analysis) | Mask_keep;
    
    curve_mean = Image(repmat(Mask_MYO,[1,1,1,nof]));
    curve_mean = reshape(curve_mean,[length(curve_mean)/nof,nof]);
    curve_mean = mean(curve_mean);
    curve_mean = curve_mean - mean(curve_mean(1:N_pre));
    pinv_curve = pinv(curve_mean)';
end





figure
plot(curve_mean)


imshow = sum(Image(:,:,:,peak_enhan_frame+(-1:1)),4);
imshow = imshow/max(imshow(:));
clear edge
for i=1:sz
    edge(:,:,i) = bwdist(Mask_MYO(:,:,i))==1;
end
imshow(edge) = 1;

ns_sqrt = ceil(sqrt(sz));
if ns_sqrt^2~=sz
    imshow(:,:,sz+1:ns_sqrt^2) = zeros(sx,sy,ns_sqrt^2-sz);
end
imshow = reshape(imshow,sx,sy*ns_sqrt,ns_sqrt);
imshow = permute(imshow,[1 3 2]);
imshow = reshape(imshow,sx*ns_sqrt,sy*ns_sqrt);
figure
imagesc(abs(imshow))
colormap gray
axis image
brighten(0.4)
drawnow
