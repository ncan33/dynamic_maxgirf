function Image_reg = registration_2D_more2more(Image,ref)
ref = abs(ref);

nof = size(Image,3);
Image_reg = Image;
parfor i=1:nof
    Image_reg(:,:,i) = Reg_GS_tv(Image(:,:,i),ref(:,:,i),1,100);
end