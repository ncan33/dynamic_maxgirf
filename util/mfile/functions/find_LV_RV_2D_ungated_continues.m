function mask = find_LV_RV_2D_ungated_continues(Image)

Image = squeeze(abs(Image));
Image_diff = diff(Image,1,3);
Image_diff_sum = sum(abs(Image_diff),3);
%figure,imagesc(Image_diff_sum)
Image_diff_sum_half = crop_half_FOV(Image_diff_sum);
Image_diff_sum_half = Image_diff_sum_half/max(Image_diff_sum_half(:));
[x,y] = find(Image_diff_sum_half==1);
x = x + size(Image_diff_sum_half,1)/2;
y = y + size(Image_diff_sum_half,2)/2;

Image_filter = zeros(size(Image_diff_sum));
Image_filter(x,y) = 1;
Image_filter = bwdist(Image_filter);
Image_filter(Image_filter==0) = 1;

mask = Image_filter<50;
%keyboard
% Image_diff_sum = Image_diff_sum./Image_filter;
% Image_diff_sum = Image_diff_sum/max(Image_diff_sum(:));

% mask = Image_diff_sum>0.05;
% mask = bwdist(mask);
% mask = mask < 5;



