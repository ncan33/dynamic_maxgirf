function Mask_all = seg_Ventricle_3D_model_based(Image,mid_slice,Mask)

[sx,sy,Nslice,nof] = size(Image);
Mask_all = false(sx,sy,Nslice);
Mask_all(:,:,mid_slice) = Mask;
order = abs((1:Nslice)-mid_slice);
[~,order]=sort(order);
order(1) = [];

curve = squeeze(sum(sum(Image(:,:,mid_slice,:).*Mask)));
Npoints = sum(sum(Mask));
curve = curve/Npoints;
[~,enhan_begin] = max(abs(diff(curve)));
N_pre_contrast = round(enhan_begin*0.8);
curve = curve - mean(curve(1:N_pre_contrast));
curve = double(curve)';

AIF_param = Quant.fit_AIF(1:nof,curve);
AIF_model = Quant.AIF_model(AIF_param,1:nof)';


[~,peak_enhan_frame] = max(smooth(curve));
Image_peak = mean(Image(:,:,mid_slice,peak_enhan_frame+(-1:1)),4);
int = mean(Image_peak(Mask));



for islice = order
    figure_temp = figure;
    flag = 1;
    flag2 = 1;
    Image_current = squeeze(Image(:,:,islice,:));
    Image_peak = Image(:,:,islice,peak_enhan_frame);
    Mask_current = false(sx,sy);
    while flag==1
        if flag2
            if islice<mid_slice
                slice_pre = islice+1;
            else
                slice_pre = islice-1;
            end
            points_add = Mask_all(:,:,slice_pre);
            flag2 = 0;
        else
            points_add =  bwdist(Mask_current)==1;
        end
        [x,y] = find(points_add);
        N = length(x);
        idx = true(1,N);
        pd_temp = zeros(1,N);
        pd_temp_int = zeros(1,N);
        pinv_AIF = pinv(AIF_model);
        for i=1:N
            curve_temp = squeeze(Image_current(x(i),y(i),:));
            curve_temp = curve_temp - mean(curve_temp(1:N_pre_contrast));
            scale = 1/(pinv_AIF*curve_temp);
            pd_temp(i) = mean(abs(AIF_model - curve_temp*scale))/max(AIF_model);
            int_temp = Image_peak(x(i),y(i));
            pd_temp_int(i) = abs((int_temp - int)/int);
        end
        idx(pd_temp>0.15|pd_temp_int>0.5) = false;
        if sum(idx)==0
            flag = 0;
        else
            x(~idx) = [];
            y(~idx) = [];
            N_add = length(x);
            Mask_add = false(sx,sy);
            for i=1:N_add
                Mask_current(x(i),y(i)) = true;
                Mask_add(x(i),y(i)) = true;
            end
            curve_temp = Image_current(repmat(Mask_add,[1,1,nof]));
            curve_temp = reshape(curve_temp,[length(curve_temp)/nof,nof]);
            Npoints_temp = size(curve_temp,1);
            curve_temp = mean(curve_temp,1)';
            curve_temp = curve_temp-mean(curve_temp(1:N_pre_contrast));
            curve = (curve*Npoints + curve_temp'*Npoints_temp)/(Npoints+Npoints_temp);
            curve = double(curve);
            AIF_param = Quant.fit_AIF(1:nof,curve);
            AIF_model = Quant.AIF_model(AIF_param,1:nof)';

            int = (sum(Image_peak(Mask_add))+int*Npoints)/(Npoints+Npoints_temp);
            Npoints = Npoints + Npoints_temp;
        end
        edge = bwdist(Mask_current)==1;
        im_show = Image_peak;
        im_show = im_show/max(im_show(:));
        im_show(edge) = 1;
        figure(figure_temp);
        clf
        imagesc(im_show)
        axis image
        colormap gray
        brighten(0.4)
        drawnow
        
    end
    close(figure_temp)
    Mask_all(:,:,islice) = imfill(Mask_current,'holes');

end




imshow = sum(Image(:,:,:,peak_enhan_frame+(-1:1)),4);
imshow = imshow/max(imshow(:));
clear edge
for i=1:Nslice
    edge(:,:,i) = bwdist(Mask_all(:,:,i))==1;
end
imshow(edge) = 1;

ns_sqrt = ceil(sqrt(Nslice));
if ns_sqrt^2~=Nslice
    imshow(:,:,Nslice+1:ns_sqrt^2) = zeros(sx,sy,ns_sqrt^2-Nslice);
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
