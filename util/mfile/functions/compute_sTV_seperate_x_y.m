function sTV_update = compute_sTV_seperate_x_y(img,weight,beta_sqrd)

if weight~=0
    diff_x = diff(img,1,2);
    diff_x = diff_x./(sqrt(beta_sqrd + abs(diff_x).^2));
    diff_x = cat(2,diff_x(:,1,:,:,:,:),diff(diff_x,1,2),-diff_x(:,end,:,:,:,:));
    
    diff_y = diff(img,1,1);
    diff_y = diff_y./(sqrt(beta_sqrd + abs(diff_y).^2));
    diff_y = cat(1,diff_y(1,:,:,:,:,:),diff(diff_y,1,1),-diff_y(end,:,:,:,:,:));

    sTV_update = weight .* (diff_x + diff_y);
else
    sTV_update = 0;
end

end