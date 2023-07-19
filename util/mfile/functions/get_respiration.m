function [phase1,phase2] = get_respiration(Image)
% [sys_loc,dia_loc] = get_sys_dia(Image)

Image = abs(Image);

figure
imagesc(sum(Image,3))
colormap gray
brighten(0.4)
axis image
axis off

h = imfreehand(gca,'closed',false); 
mask = createMask(h);

respiration_signal = squeeze(sum(sum(Image.*mask)));

phase1 = respiration_signal.'>mean(respiration_signal);
phase2 = ~phase1;