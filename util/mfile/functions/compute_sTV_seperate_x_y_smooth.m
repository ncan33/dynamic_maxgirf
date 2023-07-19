function sTV_update = compute_sTV_seperate_x_y_smooth(img,weight,beta_sqrd)

if weight~=0
    img_real = real(img);
    img_imag = imag(img);
    
    real_x = img_real(:,2:end,:,:,:,:) > img_real(:,1:end-1,:,:,:,:);
    real_x = cat(2,real_x(:,1,:,:,:,:),diff(real_x,1,2),-real_x(:,end,:,:,:,:));

    imag_x = img_imag(:,2:end,:,:,:,:) > img_imag(:,1:end-1,:,:,:,:);
    imag_x = cat(2,imag_x(:,1,:,:,:,:),diff(imag_x,1,2),-imag_x(:,end,:,:,:,:));
    
    real_y = img_real(2:end,:,:,:,:,:) > img_real(1:end-1,:,:,:,:,:);
    real_y = cat(1,real_y(1,:,:,:,:,:),diff(real_y,1,1),-real_y(end,:,:,:,:,:));

    imag_y = img_imag(2:end,:,:,:,:,:) > img_imag(1:end-1,:,:,:,:,:);
    imag_y = cat(1,imag_y(1,:,:,:,:,:),diff(imag_y,1,1),-imag_y(end,:,:,:,:,:));
    
    sTV_update = 2 .* weight .* (real_x + real_y + 1i*imag_x + 1i*imag_y);
else
    sTV_update = 0;
end

end