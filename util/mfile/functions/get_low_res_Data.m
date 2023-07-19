function Data = get_low_res_Data(Data,siz)

sx = size(Data.kSpace,1);
sy = size(Data.kSpace,2);
cut_x = (sx-siz(1))/2;
cut_y = (sy-siz(2))/2;

Data.kSpace = fftshift2(Data.kSpace);
Data.kSpace = Data.kSpace(cut_x+1:cut_x+siz(1),cut_y+1:cut_y+siz(2),:,:,:,:,:);
Data.kSpace = fftshift2(Data.kSpace);

Data.mask = fftshift2(Data.mask);
Data.mask = Data.mask(cut_x+1:cut_x+siz(1),cut_y+1:cut_y+siz(2),:,:,:,:,:);
Data.mask = fftshift2(Data.mask);

if size(Data.filter,1) > 1
    Data.filter = fftshift2(Data.filter);
    Data.filter = Data.filter(cut_x+1:cut_x+siz(1),cut_y+1:cut_y+siz(2),:,:,:,:,:);
    Data.filter = fftshift2(Data.filter);
end

Data.first_est = fft2(Data.first_est);
Data.first_est = fftshift2(Data.first_est);
Data.first_est = Data.first_est(cut_x+1:cut_x+siz(1),cut_y+1:cut_y+siz(2),:,:,:,:,:);
Data.first_est = fftshift2(Data.first_est);
Data.first_est = ifft2(Data.first_est);

Data.sens_map = fft2(Data.sens_map);
Data.sens_map = fftshift2(Data.sens_map);
Data.sens_map = Data.sens_map(cut_x+1:cut_x+siz(1),cut_y+1:cut_y+siz(2),:,:,:,:,:);
Data.sens_map = fftshift2(Data.sens_map);
Data.sens_map = ifft2(Data.sens_map);
