function mask = draw_a_mask

h = imfreehand(gca,'closed',false); 
mask=createMask(h);