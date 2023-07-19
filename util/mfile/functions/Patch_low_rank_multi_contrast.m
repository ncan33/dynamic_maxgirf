function Image = Patch_low_rank_multi_contrast(Image,patch_siz,over_lapping,idx_all)


[sx,sy,N_contrast] = size(Image);
sliding_siz = patch_siz-over_lapping;
N_patch_x = floor((sx-patch_siz(1))/sliding_siz(1));
N_patch_y = floor((sy-patch_siz(2))/sliding_siz(2));

patch_all = zeros([patch_siz,N_patch_x,N_patch_y,N_contrast]);

for i=1:N_patch_x
    for j=1:N_patch_y
        x_temp = sliding_siz(1)*(i-1)+1;
        y_temp = sliding_siz(2)*(j-1)+1;
        patch_temp = Image(x_temp:x_temp+patch_siz(1)-1,y_temp:y_temp+patch_siz(2)-1,:);
        patch_all(:,:,i,j,:) = patch_temp;
    end
end
patch_all = reshape(patch_all,[patch_siz,N_patch_x*N_patch_y,N_contrast]);


Image_low_rank = zeros(size(Image),'like',Image);
N_min = size(idx_all,3);
for i=1:N_patch_x
    i
    for j=1:N_patch_y
        idx_temp = idx_all(i,j,:);
        patch_tensor_temp = patch_all(:,:,idx_temp,:);
        patch_tensor_temp = reshape(patch_tensor_temp,[prod(patch_siz),N_min,N_contrast]);
        patch_tensor_temp = LowRank_3D_tensor(patch_tensor_temp);
        patch_tensor_temp = reshape(patch_tensor_temp(:,1,:),[patch_siz,N_contrast]);
        
        x_temp = sliding_siz(1)*(i-1)+1;
        y_temp = sliding_siz(2)*(j-1)+1;
        Image_low_rank(x_temp:x_temp+patch_siz(1)-1,y_temp:y_temp+patch_siz(2)-1,:) = Image_low_rank(x_temp:x_temp+patch_siz(1)-1,y_temp:y_temp+patch_siz(2)-1,:) + patch_tensor_temp;
    end
end
Image_low_rank(Image_low_rank==0) = Image(Image_low_rank==0);

mask = zeros(size(Image(:,:,1)),'like',Image);
for i=1:N_patch_x
    for j=1:N_patch_y
        x_temp = sliding_siz(1)*(i-1)+1;
        y_temp = sliding_siz(2)*(j-1)+1;
        mask(x_temp:x_temp+patch_siz(1)-1,y_temp:y_temp+patch_siz(2)-1) = mask(x_temp:x_temp+patch_siz(1)-1,y_temp:y_temp+patch_siz(2)-1) + 1;
    end
end
mask(mask==0) = 1;

Image = Image_low_rank./mask;