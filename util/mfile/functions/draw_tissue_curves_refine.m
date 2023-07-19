function draw_tissue_curves_refine(Image,Image_PD,mask)

Image = Image./Image_PD;
Image = Image - mean(Image(:,:,1:10),3);

N = size(mask,3);
nof = size(Image,3);
for i=1:N
    for j = 1:nof
        temp = Image(:,:,j);
        temp = temp(mask(:,:,i));
        temp_mid = median(temp);
        temp_std = std(temp);
        temp(temp<temp_mid-temp_std | temp>temp_mid+temp_std) = [];
        curve(j,i) = temp_mid;
    end
end

figure,plot(curve)
    
legend(num2str((1:6).'))