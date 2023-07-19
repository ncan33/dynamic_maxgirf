function Mask_all = seg_Ventricle_3D(Image,mid_slice,Mask)

[sx,sy,Nslice,nof] = size(Image);
Mask_all = false(sx,sy,Nslice);
Mask_all(:,:,mid_slice) = Mask;
order = abs((1:Nslice)-mid_slice);
[~,order]=sort(order);
order(1) = [];
Npoints = sum(sum(Mask));
curve = squeeze(sum(sum(Image(:,:,mid_slice,:).*Mask)))/Npoints;
%curve = curve/mean(curve);
[~,peak_enhan_frame] = max(smooth(curve));
%Image_peak = mean(Image(:,:,mid_slice,peak_enhan_frame+(-1:1)),4);
%int = mean(Image_peak(Mask));

pinv_curve = pinv(curve);

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
        for i=1:N
            curve_temp = squeeze(Image_current(x(i),y(i),:));
            scale = 1/(pinv_curve*curve_temp);
            %curve_temp = curve_temp/mean(curve_temp);
            pd_temp = mean(abs(curve - curve_temp*scale)./curve);
            %int_temp = Image_peak(x(i),y(i));
            %pd_temp_int = abs((int_temp - int)/int);
            if pd_temp > 0.1 %|| pd_temp_int >0.5
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
                Mask_current(x(i),y(i)) = true;
            end
            curve_temp = Image_current(repmat(Mask_current,[1,1,nof]));
            curve_temp = reshape(curve_temp,[length(curve_temp)/nof,nof]);
            Npoints_temp = size(curve_temp,1);
            curve_temp = mean(curve_temp)';
            %curve_temp = curve_temp/mean(curve_temp);
            curve = (curve*Npoints + curve_temp*Npoints_temp)/(Npoints+Npoints_temp);
            pinv_curve = pinv(curve);
            %curve = curve/mean(curve);

            %int = (sum(Image_peak(Mask_current))+int*Npoints)/(Npoints+Npoints_temp);
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
    
    Mask_all(:,:,islice) = Mask_current;
end