function draw_ROI_binary(Image,mask,colors)
figure(1)
clf
Image = Image/max(Image(:));
Image = repmat(Image,[1,1,3]);
for seg = 1:size(mask,3)
    mask_temp = mask(:,:,seg);
    Image = Image.*mask_temp.*permute(colors(seg,:),[3,1,2])*1.3 + Image.*~mask_temp;
end
imagesc(Image),axis image, colormap gray,

for seg = 1:size(mask,3)
    patch([0,0,0,0],[0,0,0,0],colors(seg,:),'EdgeColor','none')
end

legend(num2str((1:size(mask,3)).'))
axis([1,size(Image,2),1,size(Image,1)])
axis off
end