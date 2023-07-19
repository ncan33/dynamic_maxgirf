function img = LLR_2D_STCR_ADMM(Data,para)

disp('Performing ADMM reconstruction...');
disp('Showing progress...')

Data.sens_map_conj = conj(Data.sens_map);

weight_ttv = para.Recon.weight_tTV;

para.Recon.weight_tTV = 0;
para.Recon.weight_l2 = 0.1;
para.Recon.noi = 20;
img = STCR_conjugate_gradient_ADMM(Data,para);

u = zeros(size(img), 'single');

for i=1:10
    z = sfth(circshift(img,[0,0,1]) + u, weight_ttv/0.1);
    u = u + circshift(img,[0,0,1]) - z;
    Data.first_est = img;
    Data.first_guess = u;
    img = STCR_conjugate_gradient_ADMM(Data,para);
    
    
end


end
