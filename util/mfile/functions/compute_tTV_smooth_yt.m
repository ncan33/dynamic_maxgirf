function tTV_update = compute_tTV_smooth_yt(image,weight,beta_square)
%tTV_update = compute_tTV_yt(image,weight,beta_square)
if weight~=0
    image_real = real(image);
    image_imag = imag(image);

    temp_b_real = image_real(:,:,2:end,:,:) > image_real(:,:,1:end-1,:,:);
    temp_b_imag = image_imag(:,:,2:end,:,:) > image_imag(:,:,1:end-1,:,:);
    
    temp_c_real = diff(temp_b_real,1,3);
    temp_c_imag = diff(temp_b_imag,1,3);
   
    tTV_update_real = weight * cat(3,temp_b_real(:,:,1,:,:,:),temp_c_real,-temp_b_real(:,:,end,:,:,:));
    tTV_update_imag = weight * cat(3,temp_b_imag(:,:,1,:,:,:),temp_c_imag,-temp_b_imag(:,:,end,:,:,:));
    tTV_update = tTV_update_real + 1i*tTV_update_imag;
    
    tTV_update = 2*tTV_update;
else
    tTV_update = 0;
end

end