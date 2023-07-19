function plot_curve(Image)

figure
imagesc(sum(Image,3))
colormap gray
brighten(0.4)
axis image
axis off

h = imfreehand(gca,'closed',false); 
mask=createMask(h);
close

value = squeeze(sum(sum(mask.*Image,1),2));
figure,plot(value/max(value))
