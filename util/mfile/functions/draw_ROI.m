function draw_ROI(Image,mask,colors)
figure(1)
clf
imagesc(Image),axis image, colormap gray, brighten(0.4),
hold on
for seg = 1:size(mask,3)
    [y,x] = find(mask(:,:,seg));
    plot(x,y,'.','Color',colors(seg,:),'MarkerSize',15)
end
legend(num2str((1:size(mask,3)).'))
end