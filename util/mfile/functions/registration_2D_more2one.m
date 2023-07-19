function Im_reg = registration_2D_more2one(Image,ref)
ref = abs(ref);
Image_reg = imgaussfilt(abs(Image),1);
nof = size(Image,3);
for i=1:nof
    [~,x,y,~]= Reg_GS_accurate(Image_reg(:,:,i),ref,0.1,100);
    Im_reg(:,:,i) = interp2(Image(:,:,i),y,x);
end