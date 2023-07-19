function image = compute_l1_all(image,weight)
%tTV_update = compute_tTV_yt(image,weight,beta_square)
if weight~=0
    image_r = real(image);
    image_i = imag(image);
    image_r = mean(image_r,3) - image_r;
    image_i = mean(image_i,3) - image_i;
    image_r = sign(image_r);
    image_i = sign(image_i);
    image = image_r + 1i*image_i;
    image = weight*image;
else
    image = 0;
end

end