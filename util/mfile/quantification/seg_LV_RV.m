function [mask_LV,mask_RV] = seg_LV_RV(Image,LV,RV)


Image = abs(squeeze(Image));
[sx,sy,nof] = size(Image);

LV_curve = squeeze(sum(sum(Image(LV(2)+(-1:1),LV(1)+(-1:1),:))))/9;
%LV_curve = LV_curve/mean(LV_curve);

mask_LV = false(sx,sy);
mask_LV(LV(2)+(-1:1),LV(1)+(-1:1)) = true;

[~,peak_enhan_frame] = max(smooth(LV_curve));
Image_peak = mean(Image(:,:,peak_enhan_frame + (-1:1)),3);

int_LV = mean(Image_peak(mask_LV));
figure_LV = figure;

flag = 1;
while flag == 1
    points_add = bwdist(mask_LV)==1;
    [x,y] = find(points_add);
    N = length(x);
    idx = true(1,N);
    for i=1:N
        curve_temp = squeeze(Image(x(i),y(i),:));
        %curve_temp = curve_temp/mean(curve_temp);
        pd_temp = mean(abs(LV_curve - curve_temp)./LV_curve);
        int_temp = Image_peak(x(i),y(i));
        pd_temp_int = abs((int_temp - int_LV)/int_LV);
        
        if pd_temp > 0.1 % || pd_temp_int >0.5
            idx(i) = false;
        end
        
    end
    if sum(idx)==0
        flag = 0;
    else
        x(~idx) = [];
        y(~idx) = [];
        N_add = length(x);
        for i=1:N_add
            mask_LV(x(i),y(i)) = true;
        end
        LV_curve = Image(repmat(mask_LV,[1,1,nof]));
        LV_curve = reshape(LV_curve,[length(LV_curve)/nof,nof]);
        LV_curve = mean(LV_curve,1)';
        %LV_curve = LV_curve/mean(LV_curve);
        int_LV = mean(Image_peak(mask_LV));
    end
    edge = bwdist(mask_LV)==1;
    im_show = Image_peak;
    im_show = im_show/max(im_show(:));
    im_show(edge) = 1;
    figure(figure_LV);
    clf
    imagesc(im_show)
    axis image
    colormap gray
    brighten(0.4)
    drawnow
end

return


RV_curve = squeeze(sum(sum(Image(RV(2)+(-1:1),RV(1)+(-1:1),:))));
RV_curve = RV_curve/mean(RV_curve);

mask_RV = false(sx,sy);
mask_RV(RV(2)+(-1:1),RV(1)+(-1:1)) = true;

[~,peak_enhan_frame] = max(smooth(RV_curve));
Image_peak = mean(Image(:,:,peak_enhan_frame + (-1:1)),3);

int_RV = mean(Image_peak(mask_RV));
figure_RV = figure;


flag = 1;
while flag == 1
    points_add = bwdist(mask_RV)==1;
    [x,y] = find(points_add);
    N = length(x);
    idx = true(1,N);
    for i=1:N
        curve_temp = squeeze(Image(x(i),y(i),:));
        curve_temp = curve_temp/mean(curve_temp);
        pd_temp = mean(abs(RV_curve - curve_temp)./RV_curve);
        int_temp = Image_peak(x(i),y(i));
        pd_temp_int = abs((int_temp - int_RV)/int_RV);
        
        if pd_temp > 0.2 || pd_temp_int >0.5
            idx(i) = false;
        end
        
    end
    if sum(idx)==0
        flag = 0;
    else
        x(~idx) = [];
        y(~idx) = [];
        N_add = length(x);
        for i=1:N_add
            mask_RV(x(i),y(i)) = true;
        end
        RV_curve = Image(repmat(mask_RV,[1,1,nof]));
        RV_curve = reshape(RV_curve,[length(RV_curve)/nof,nof]);
        RV_curve = mean(RV_curve,1)';
        RV_curve = RV_curve/mean(RV_curve);
        int_RV = mean(Image_peak(mask_RV));
    end
    edge = bwdist(mask_RV)==1;
    im_show = Image_peak;
    im_show = im_show/max(im_show(:));
    im_show(edge) = 1;
    figure(figure_RV);
    clf
    imagesc(im_show)
    axis image
    colormap gray
    brighten(0.4)
    drawnow
end
