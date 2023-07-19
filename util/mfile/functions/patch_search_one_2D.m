function offset = patch_search_one_2D(patch,search_window)

patch_size = size(patch);
search_size = size(search_window);
Npatches = search_size - patch_size + 1;
% Npatches_all = prod(Npatches);

patches_all = zeros([patch_size,Npatches]);

for i=1:Npatches(1)
    for j=1:Npatches(2)
            patches_all(:,:,i,j) = search_window(i:i+patch_size(1)-1,j:j+patch_size(2)-1);
    end
end

d = patch - patches_all;
d = sum(sum(abs(d)));
[~,idx] = min(d(:));
[x,y] = ind2sub(Npatches,idx);
offset = [x,y];

end