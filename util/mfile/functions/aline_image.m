function [Image_all,seed_all] = aline_image(Image_all,seed_all)

nof = size(Image_all,3);

center_all_x = seed_all(1,1,:);
center_all_y = seed_all(2,1,:);

center_mean_x = round(mean(center_all_x));
center_mean_y = round(mean(center_all_y));

dx_all = center_mean_x - center_all_x;
dy_all = center_mean_y - center_all_y;

for i=1:nof
    Image_all(:,:,i) = imtranslate(Image_all(:,:,i),[dx_all(i),dy_all(i)]);
end

seed_all = seed_all + [dx_all;dy_all];
