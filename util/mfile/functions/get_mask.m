function mask = get_mask(Image)
%Image = abs(Image);

temp = figure;
set(temp,'Position',[10,10,1000,1000]);
imagesc(sum(squeeze(Image),3));
colormap gray
brighten(0.4)
axis image
axis off
h = imfreehand(gca,'closed',false); 
mask = createMask(h);
close(temp)