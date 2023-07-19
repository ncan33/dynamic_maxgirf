function show_patch(image,patch_loc,patch_size,offset,slice)


patch_s = patch_loc;
patch_e = patch_loc + patch_size - 1;
color = colors;
figure
for i=1:size(image,5)
    
    subplot(1,2,1)
    
    patch_s = patch_s + offset(i,:);
    patch_e = patch_e + offset(i,:);
    
    im = squeeze(abs(image(:,:,patch_s(3)+slice(1)-1,patch_s(4)+slice(2)-1,i)));
    imagesc(im)
    colormap gray
    axis image
    brighten(0.4)
    title(i)
    hold on
    drawnow

    
    plot([patch_s(1),patch_e(1),patch_e(1),patch_s(1),patch_s(1)],[patch_e(2),patch_e(2),patch_s(2),patch_s(2),patch_e(2)],'Color',color(1,:))
    
    subplot(1,2,2)
    patch = image(patch_s(1):patch_s(1)+patch_size(1)-1,patch_s(2):patch_s(2)+patch_size(2)-1,patch_s(3)+slice(1)-1,patch_s(4)+slice(2)-1,i);
    imagesc(patch)
    colormap gray
    axis image
    brighten(0.4)
    drawnow

end