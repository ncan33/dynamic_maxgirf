function Image = LLR_2D_ADMM(Data,para)

disp('Performing ADMM reconstruction...');
disp('Showing progress...')


para.setting.ifGPU = 0;
para.Recon.noi = 30;

para.Recon.weight_l2 = 0.1;

Image = STCR_conjugate_gradient_ADMM(Data,para);
Data.Y = zeros(size(Image));

tau = 3000;

for i = 1:5
    
    Image_lr = low_rank_yt(Image + Data.Y, 8, 8, tau);
    
    Data.Y = Data.Y + Image - Image_lr;
    
    Data.first_guess = Image_lr;
    Data.first_est = Image;

    Image = STCR_conjugate_gradient_ADMM(Data,para);
end


end