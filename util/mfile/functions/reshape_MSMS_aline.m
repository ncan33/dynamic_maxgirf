function Image = reshape_MSMS_aline(Image)
siz = size(Image);
siz5 = size(Image,5);
Image = permute(Image,[1,2,3,5,4]);
Image = Image(:,:,:,:);
Image = permute(Image,[1,2,4,3]);
Image = reshape(Image,siz(1),siz(2)*siz(4)*siz5,siz(3));