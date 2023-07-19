addpath /v/raid1b/ytian/else/ModelBasedMy/

nSlice = size(Image,3);
Image_reg = zeros(size(Image));

for i=1:nSlice
    slice_temp = squeeze(Image(:,:,i,:));
    %PD_slcie_temp = Image_PD(:,:,i);
    slice_temp = rigid_regisgtration_2D(slice_temp);
    MBI_temp = run_registration(slice_temp);
    slice_temp = registration_2D_more2more(slice_temp,MBI_temp);
    Image_reg(:,:,i,:) = slice_temp;
end


Image_normalized = Image_reg./Image_PD;
Image_normalized = Image_normalized(:,:,:,6:end);
Image_normalized = Image_normalized - mean(Image_normalized(:,:,:,1:10),4);

