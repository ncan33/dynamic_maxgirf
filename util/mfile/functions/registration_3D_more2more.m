function Im_reg = registration_3D_more2more(Image,ref)
ref = abs(ref);

nof = size(Image,4);
for i=1:nof
    tic
    Im_reg(:,:,:,i) = Reg_GS_3D_new(Image(:,:,:,i),ref(:,:,:,i),1,20);
    toc
end