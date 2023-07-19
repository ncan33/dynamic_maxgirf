function tTV_update = compute_3DtTV_circ(image,weight,beta_square)

if weight~=0
    image = cat(4, image(:, :, :, end, :, :), image);
    tTV_update = diff(image,1,4);
    tTV_update = tTV_update./(sqrt(beta_square+(abs(tTV_update).^2)));
    tTV_update = cat(4, tTV_update(:, :, :, end, :, :), tTV_update);
    tTV_update = diff(tTV_update,1,4);
    tTV_update = weight * tTV_update;
    tTV_update = circshift(tTV_update, -1, 4);
else
    tTV_update = 0;
end

end