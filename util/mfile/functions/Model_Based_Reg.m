function Image_reg = Model_Based_Reg(Image)

% addpath /v/raid1b/ytian/else/ModelBasedMy/

noi = 5;
Image_reg = Image;
for i=1:noi
    Image_MB = run_registration(Image_reg);
    Image_reg = registration_2D_more2more(Image,Image_MB);
end
