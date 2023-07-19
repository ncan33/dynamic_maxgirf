function Im_reg = registration_2D_edge_one2more(Image,ref)
ref_edge = abs(get_total_variation(ref));
Image_edge = abs(get_total_variation(Image));
nof = size(Image,3);
for i=1:nof
    [~,x,y,~]= Reg_GS_accurate(ref_edge,Image_edge(:,:,i),0.1,100);
    Im_reg(:,:,i) = interp2(ref,y,x);
end