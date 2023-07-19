function target = rescale_for_PSF(source,target,mask)

[~,~,nof,nSMS,ns] = size(source);

for i=1:nof
    for j=1:nSMS
        for k=1:ns
            temp_s = source(:,:,i,j,k);
            temp_t = target(:,:,i,j,k);
            temp_m = mask(:,:,i,j,k);
            temp_scale = abs(temp_s(:)).' * pinv(abs(temp_t(:).*temp_m(:))).';
            target(:,:,i,j,k) = target(:,:,i,j,k)*temp_scale;
        end
    end
end

            