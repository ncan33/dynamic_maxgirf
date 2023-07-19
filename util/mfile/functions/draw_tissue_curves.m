function draw_tissue_curves(Image,Image_PD,mask,TimeStamp_second)

Image = Image./Image_PD;
Image = Image - mean(Image(:,:,1:10),3);

N = size(mask,3);
for i=1:N
    curve(:,i) = squeeze(sum(sum(Image.*mask(:,:,i))))/sum(sum(mask(:,:,i)));
end

figure,plot(TimeStamp_second-TimeStamp_second(1),curve,'LineWidth',1.5)
    
legend(num2str((1:N).'))