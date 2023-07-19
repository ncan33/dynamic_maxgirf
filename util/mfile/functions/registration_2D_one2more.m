function Im_reg = registration_2D_one2more(Image,ref)
ref = abs(ref);
Image_reg = imgaussfilt(abs(Image),1);
nof = size(Image,3);
for i=1:nof
    [~,x,y,~]= Reg_GS_accurate(ref,Image_reg(:,:,i),0.1,100);
    Im_reg(:,:,i) = interp2(ref,y,x);
end