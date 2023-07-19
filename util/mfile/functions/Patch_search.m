function idx_all = Patch_search(Image,patch_siz,over_lapping)

[sx,sy] = size(Image);
sliding_siz = patch_siz-over_lapping;
N_patch_x = floor((sx-patch_siz(1))/sliding_siz(1));
N_patch_y = floor((sy-patch_siz(2))/sliding_siz(2));

patch_all = zeros([patch_siz,N_patch_x,N_patch_y]);

for i=1:N_patch_x
    for j=1:N_patch_y
        x_temp = sliding_siz(1)*(i-1)+1;
        y_temp = sliding_siz(2)*(j-1)+1;
        patch_temp = Image(x_temp:x_temp+patch_siz(1)-1,y_temp:y_temp+patch_siz(2)-1);
        patch_all(:,:,i,j) = patch_temp;
    end
end



N_min = 30;
idx_all = zeros(N_patch_x,N_patch_y,N_min);
for i=1:N_patch_x
    for j=1:N_patch_y
        d_temp = patch_all(:,:,i,j) - patch_all;
        d_temp = abs(d_temp);
        d_temp = squeeze(sum(sum(d_temp)));
        d_temp = d_temp(:);
        [~,idx_temp] = mink(d_temp,N_min);
        idx_all(i,j,:) = idx_temp;
    end
end


%patch_all = reshape(patch_all,[patch_siz,N_patch_x*N_patch_y]);