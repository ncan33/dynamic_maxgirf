function Mask = seg_Ventricle_3D_step_2(Image,Mask_all)

[sx,sy,sz,nof] = size(Image);
curves = Image(repmat(Mask_all,[1,1,1,nof]));
Npoints = length(curves)/nof;

curves = reshape(curves,[Npoints,nof]);
curve_mean = mean(curves)';
[~,peak_enhan_frame] = max(smooth(curve_mean));
[~,enhan_begin] = max(abs(diff(smooth(curve_mean))));
N_pre = round(enhan_begin*0.8);

% curve_mean = curve_mean-mean(curve_mean(1:N_pre));
% AIF_param = Quant.fit_AIF(1:nof,double(curve_mean));
% AIF_model = Quant.AIF_model(AIF_param,1:nof)';

Mask = Mask_all;

Niter = 10;
for iter=1:Niter
    Mask_analysis = bwdist(Mask) == 1 | bwdist(~Mask) == 1;
    
    idx = find(Mask_analysis);
    [x,y,z] = ind2sub([sx,sy,sz],idx);
    N_analysis = length(x);
    %pinv_AIF = pinv(AIF_model);
    pinv_curve = pinv(curve_mean);
    pd_temp = zeros(N_analysis,1);
    for i=1:N_analysis
        curve_temp = squeeze(Image(x(i),y(i),z(i),:));
        %curve_temp = curve_temp - mean(curve_temp(1:N_pre));
        scale = 1/(pinv_curve*curve_temp);
        pd_temp(i) = mean(abs(curve_mean - curve_temp*scale)./curve_mean);
    end
    keep = pd_temp<0.1;
    
    Mask_keep = false(sx,sy,sz);
    Mask_keep(idx(keep)) = true;
    Mask = Mask_all & ~Mask_analysis | Mask_keep;
    
    for i=1:sz
        Mask(:,:,i) = bwareafilt(Mask(:,:,i),1);
        Mask(:,:,i) = bwfill(Mask(:,:,i),'holes');
    end    
    
    curves = Image(repmat(Mask,[1,1,1,nof]));
    Npoints = length(curves)/nof;
    curves = reshape(curves,[Npoints,nof]);
    curve_mean = mean(curves)';
%     curve_mean = curve_mean-mean(curve_mean(1:N_pre));
%     AIF_param = Quant.fit_AIF(1:nof,double(curve_mean));
%     AIF_model = Quant.AIF_model(AIF_param,1:nof)';
    
end

imshow = sum(Image(:,:,:,peak_enhan_frame+(-1:1)),4);
imshow = imshow/max(imshow(:));
clear edge
for i=1:sz
    edge(:,:,i) = bwdist(Mask(:,:,i))==1;
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
