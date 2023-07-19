function tTV_update = compute_tTV_random(image,weight,beta_square)
%tTV_update = compute_tTV_yt(image,weight,beta_square)

nof = size(image,3);
order = rand(1,nof);
[~,order] = sort(order);
[~,order_back] = sort(order);
image = image(:,:,order);

if weight~=0
    image = image(:,:,[end,1:end,1]);
    image = diff(image,1,3);
    image = image./(sqrt(beta_square+(abs(image).^2)));
    tTV_update = weight .* diff(image,1,3);
%     tTV_update =  cat(3,temp_b(:,:,1,:,:,:),temp_c,-temp_b(:,:,end,:,:,:));
else
    tTV_update = 0;
end

tTV_update = tTV_update(:,:,order_back);

end