function slTV_update = compute_slice_TV_MSMS(image,weight,beta_square)

if weight~=0
    siz = size(image);
    nslice = siz(5)*siz(6);
    order = 1:siz(5):nslice;
    for i=siz(5):-1:2
        order = [order,i:siz(5):nslice];
    end
    [~,order_back] = sort(order);
    image = image(:,:,:,:,order);
    temp_a = diff(image,1,5);
    temp_b = temp_a./(sqrt(beta_square+(abs(temp_a).^2)));
    temp_c = diff(temp_b,1,5);
    slTV_update = weight .* cat(5,temp_b(:,:,:,:,1,:),temp_c,-temp_b(:,:,:,:,end,:));
    slTV_update = slTV_update(:,:,:,:,order_back);
    slTV_update = reshape(slTV_update,siz);
else
    slTV_update = zeros(size(image),'like',image);
end

end