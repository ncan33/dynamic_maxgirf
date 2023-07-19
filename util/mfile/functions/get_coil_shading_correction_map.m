function CSCMap = get_coil_shading_correction_map(Image)

sx = size(Image,1);
sy = size(Image,2);
ns = size(Image,3);
Image = abs(Image);
CSCMap = zeros(size(Image),'like',Image);


for i=1:ns
    Image_temp = double(Image(:,:,i));
    regional_max = imregionalmax(Image_temp);
    MaxInt = max(Image_temp(:));
    Threshold = 0.1*MaxInt;
    %Image_temp(Image_temp<Threshold) = Threshold;
    mask = Image_temp > 0;
    mask = mask.*regional_max;
    %mask_edge = false(size(mask));
    %mask_edge(1,1:end) = true;
    %mask_edge(end,1:end) = true;
    %mask_edge(1:end,1) = true;
    %mask_edge(1:end,end) = true;
    %mask_edge = mask_edge.*~mask;
    idx = find(mask);
    [x,y] = ind2sub(size(Image_temp),idx);keyboard
    [Y,X] = meshgrid(1:sx,1:sy);
    %map_nearest = griddata(x,y,Image_temp(idx),X,Y,'nearest');
    map_temp = griddata(x,y,Image_temp(idx),X,Y,'cubic');
    CSCMap(:,:,i) = map_temp;
end