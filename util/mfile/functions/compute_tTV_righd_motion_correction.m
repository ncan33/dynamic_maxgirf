function tTV_update = compute_tTV_righd_motion_correction(image,Motion,weight,beta_square)

if weight~=0
    temp_a = diff(image,1,3);
    temp_b = temp_a./(sqrt(beta_square+(abs(temp_a).^2)));
    temp_c = diff(temp_b,1,3);
    tTV_update = weight * cat(3,temp_b(:,:,1,:,:,:),temp_c,-temp_b(:,:,end,:,:,:));
    
    
    ROI_update = zeros(size(tTV_update),'like',tTV_update);
    for i = Motion.slice_pick
        temp = image(:,:,:,:,i);
        temp_shift = Motion.y_motion(:,i);
        for t = 1:size(ROI_update,3)
            temp(:,:,t) = circshift(temp(:,:,t),-temp_shift(t),1);
        end
        temp_a = diff(temp,1,3);
        temp_b = temp_a./(sqrt(beta_square+(abs(temp_a).^2)));
        temp_c = diff(temp_b,1,3);
        temp = weight * cat(3,temp_b(:,:,1,:,:,:),temp_c,-temp_b(:,:,end,:,:,:));
        for t=1:size(ROI_update,3)
            temp(:,:,t) = circshift(temp(:,:,t),temp_shift(t),1);
        end
        ROI_update(:,:,:,:,i) = temp;
    end
    
    tTV_update(Motion.Mask) = ROI_update(Motion.Mask);

else
    tTV_update = 0;
end

end