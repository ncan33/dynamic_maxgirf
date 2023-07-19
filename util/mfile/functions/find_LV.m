function [LV,RV] = find_LV(Image)
% [LV,RV] = find_LV(Image)

Image = abs(squeeze(Image));
Image_small = (Image);
Image_std = std(Image_small,1,3);
Image_max = max(Image_small,[],3);
Image_std = imgaussfilt(Image_std,3);
Image_max = imgaussfilt(Image_max,3);

Image_std = Image_max.*Image_std;

image_center = size(Image_std)/2;
image_center = image_center(1:2);


threshold = 0.5+0.01;

x = [];
y = [];

while length(x)<2
    threshold = threshold - 0.01;
    Mask = imregionalmax(Image_std) & Image_std>max(Image_std(:)*threshold);
    [x,y] = find(Mask);
    distance_to_center = sqrt((y-image_center(1)).^2+(x-image_center(2)).^2);
    idx_drop = distance_to_center > mean(image_center)*0.45;
    %distance_to_center(idx_drop) = [];
    x(idx_drop) = [];
    y(idx_drop) = [];
    curves = zeros(length(x),size(Image,3));
    for i=1:length(x)
        curves(i,:) = Image(x(i),y(i),:);
    end
    max_curves = max(curves,[],2);
    idx_low_int = max_curves<max(max_curves)*0.5;
    %curves(idx_low_int,:) = [];
    x(idx_low_int) = [];
    y(idx_low_int) = [];
    
end

if length(x)==2
    curve1 = sum(sum(Image_small(x(1)-1:x(1)+1,y(1)-1:y(1)+1,:)));
    curve2 = sum(sum(Image_small(x(2)-1:x(2)+1,y(2)-1:y(2)+1,:)));
    [peak1_value,peak1] = max(curve1);
    [peak2_value,peak2] = max(curve2);
    if peak1_value*0.6 > peak2_value
        LV = [y(1),x(1)];
        RV = [y(2),x(2)];
        T = 'LV:red   RV:black';
    elseif peak2_value*0.6 > peak1_value
        LV = [y(2),x(2)];
        RV = [y(1),x(1)];
        T = 'LV:red   RV:black';
    elseif peak1>size(Image,3)/2*0.9
        RV = [y(1),x(1)];
        LV = [y(2),x(2)];
        T = 'LV:red   RV:black';
    elseif peak2>size(Image,3)/2*0.9
        RV = [y(2),x(2)];
        LV = [y(1),x(1)];
        T = 'LV:red   RV:black';
        
    else
        
        if peak1<peak2
            LV = [y(2),x(2)];
            RV = [y(1),x(1)];
            T = 'LV:red   RV:black';
        elseif peak1==peak2
            LV = [y(1),x(1)];
            RV = [y(2),x(2)];
            warning 'LV/RV might be both LV/RV'
            T = 'LV/RV might be both LV/RV';
        else
            LV = [y(1),x(1)];
            RV = [y(2),x(2)];
            T = 'LV:red   RV:black';
        end
        figure
    end
    imagesc(Image_max);hold on; axis image;
    title(T)
    plot(LV(1),LV(2),'ro','MarkerFaceColor','r')
    plot(RV(1),RV(2),'ko','MarkerFaceColor','k')

else
    curves = zeros(length(x),size(Image,3));
    for i=1:length(x)
        curves(i,:) = sum(sum(Image(x(i)+(-1:1),y(i)+(-1:1),:)));
    end
    max_curves = max(curves,[],2);
    idx_low_int = max_curves<max(max_curves)*0.5;
    curves(idx_low_int,:) = [];
    x(idx_low_int) = [];
    y(idx_low_int) = [];
    bins = kmeans(curves,2);
    %curve = zeros(2,size(Image,3));
    idx1 = bins==1;
    x1 = x(~idx1);   
    y1 = y(~idx1);
    curve1 = curves(~idx1,:);
    x2 = x(idx1);
    y2 = y(idx1);
    curve2 = curves(idx1,:);
    
    [~,curve1_max] = max(max(curve1,[],2));
    x1 = x1(curve1_max);
    y1 = y1(curve1_max);
    [~,curve2_max] = max(max(curve2,[],2));
    x2 = x2(curve2_max);
    y2 = y2(curve2_max);
    
    curve1_mean = mean(curve1,1);
    curve2_mean = mean(curve2,1);

    [peak1_value,peak1] = max(curve1_mean);
    [peak2_value,peak2] = max(curve2_mean);
    if peak1_value*0.6 > peak2_value
        LV = [y1,x1];
        RV = [y2,x2];
        T = 'LV:red   RV:black';
    elseif peak2_value*0.6 > peak1_value
        LV = [y2,x2];
        RV = [y1,x1];
        T = 'LV:red   RV:black';
    elseif peak1>size(Image,3)/2*0.9
        RV = [y1,x1];
        LV = [y2,x2];
        T = 'LV:red   RV:black';
    elseif peak2>size(Image,3)/2*0.9
        RV = [y2,x2];
        LV = [y1,x1];
        T = 'LV:red   RV:black';
    else
        if peak1<peak2
            LV = [y2,x2];
            RV = [y1,x1];
            T = 'LV:red   RV:black';
        elseif peak1==peak2
            if peak1_value>peak2_value
                LV = [y1,x1];
                RV = [y2,x2];
            else
                RV = [y1,x1];
                LV = [y2,x2];
            end
            warning 'LV/RV might be both LV/RV'
            T = 'LV/RV might be both LV/RV';
        else
            LV = [y1,x1];
            RV = [y2,x2];
            T = 'LV:red   RV:black';
        end
    end
    figure
    imagesc(Image_max);hold on; axis image;
    title(T)
    plot(LV(1),LV(2),'ro','MarkerFaceColor','r')
    plot(RV(1),RV(2),'ko','MarkerFaceColor','k')
        
end
