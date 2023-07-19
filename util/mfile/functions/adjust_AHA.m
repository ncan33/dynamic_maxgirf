function B = adjust_AHA(B)

clc

Nslice = size(B.bld_flow_map,3);


for i=1:Nslice
    close all
    background_image = B.Image_show(:,:,i);
    color_image = B.bld_flow_map(:,:,i).*B.Mask_all(:,:,i);
    background_image = background_image/max(background_image(:));
    background_image = repmat(background_image,[1,1,3]);
    
    temp = figure;
    set(temp,'Position',[10,10,1000,1000]);
    imagesc(color_image)
    colormap jet
    hold on
    bg = imagesc(background_image);
    
    mask = logical(color_image);
    set(bg,'AlphaData',~mask);
    
    colorbar
    colormap hot
    axis image
    axis off
    
    title(['Slice ',num2str(i)])
    
    prompt = 'Do I want to adjust the mask? (0/1) ';
    adjust_flag = input(prompt);
    
    while adjust_flag
        h = imfreehand(gca,'closed',false);
        mask_delete = createMask(h,bg);
        B.Seg_mask{i} = B.Seg_mask{i}.*~mask_delete;
        B.Mask_all(:,:,i) = B.Mask_all(:,:,i).*~mask_delete;
        close all
        background_image = B.Image_show(:,:,i);
        color_image = B.bld_flow_map(:,:,i).*B.Mask_all(:,:,i);
        background_image = background_image/max(background_image(:));
        background_image = repmat(background_image,[1,1,3]);
        
        temp = figure;
        set(temp,'Position',[10,10,1000,1000]);
        imagesc(color_image)
        colormap jet
        hold on
        bg = imagesc(background_image);
        
        mask = logical(color_image);
        set(bg,'AlphaData',~mask);
        
        colorbar
        colormap hot
        axis image
        axis off
        
        title(['Slice ',num2str(i)])
        prompt = 'Do I want to adjust the mask? (0/1) ';
        adjust_flag = input(prompt);
    end
    
end
keyboard

for i=1:Nslice
    B.Seg_mask{i} = logical(B.Seg_mask{i});
    Nseg = size(B.Seg_mask{i},3);
    bld_flow_map_temp = B.bld_flow_map(:,:,i);
    for j=1:Nseg
        mask_temp = B.Seg_mask{i}(:,:,j);
        bld_flow_temp = bld_flow_map_temp(mask_temp);
        bld_flow(i,j,1) = mean(bld_flow_temp);
        bld_flow(i,j,2) = std(bld_flow_temp);
    end
end

B.bld_flow = bld_flow;

close all
for i=1:Nslice
    
    background_image = B.Image_show(:,:,i);
    color_image = B.bld_flow_map(:,:,i).*B.Mask_all(:,:,i);
    background_image = background_image/max(background_image(:));
    background_image = repmat(background_image,[1,1,3]);
    
    temp = figure;
    set(temp,'Position',[10,10,1000,1000]);
    imagesc(color_image)
    colormap jet
    hold on
    bg = imagesc(background_image);
    
    mask = logical(color_image);
    set(bg,'AlphaData',~mask);
    
    colorbar
    colormap hot
    axis image
    axis off
    
    title(['Slice ',num2str(i)])
    
end
% B.Seg_mask = Seg_mask;
% B.Mask_all = Mask_all;
% B.Image_show = Image_peak;
% B.bld_flow_map = bld_flow_map;