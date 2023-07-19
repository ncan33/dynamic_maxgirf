function Ref = get_4_reference(Image)

Image_temp = Image.bin1;
ref = squeeze(mean(Image_temp,3));
for i=1:size(ref,3)
    reg_temp(:,:,:,i) = registration_2D(Image_temp(:,:,:,i),ref(:,:,i));
end
Ref.bin1 = abs(mean(reg_temp,3));clear reg_temp;

Image_temp = Image.bin2;
ref = squeeze(mean(Image_temp,3));
for i=1:size(ref,3)
    reg_temp(:,:,:,i) = registration_2D(Image_temp(:,:,:,i),ref(:,:,i));
end
Ref.bin2 = abs(mean(reg_temp,3));clear reg_temp;

Image_temp = Image.bin3;
ref = squeeze(mean(Image_temp,3));
for i=1:size(ref,3)
    reg_temp(:,:,:,i) = registration_2D(Image_temp(:,:,:,i),ref(:,:,i));
end
Ref.bin3 = abs(mean(reg_temp,3));clear reg_temp;

Image_temp = Image.bin4;
ref = squeeze(mean(Image_temp,3));
for i=1:size(ref,3)
    reg_temp(:,:,:,i) = registration_2D(Image_temp(:,:,:,i),ref(:,:,i));
end
Ref.bin4 = abs(mean(reg_temp,3));clear reg_temp;