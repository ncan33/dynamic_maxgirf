function tTV_update = compute_tTV_bins(Image,weight,beta_square,bins)
tTV_update = zeros(size(Image),'like',Image);

for i=1:size(bins,1)
    bin_temp = bins(i,:);
    Image_temp = Image(:,:,bin_temp,:,:);
    
    tTV_update(:,:,bin_temp,:,:) = compute_tTV_yt(Image_temp,weight,beta_square);
    tTV_update(:,:,bin_temp,:,:) = tTV_update(:,:,bin_temp,:,:) + (LowRank_yt(Image_temp)-Image_temp)*0.05;
end