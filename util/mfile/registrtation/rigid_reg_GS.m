function [Image_reg,shifts,mask] = rigid_reg_GS(Image,mask)

if ~exist('mask')
    mask = get_mask(Image(:,:,end));
end
Image_reg = Image;
nof = size(Image,3);
noi = 4;
step = 100;

shifts = zeros(2,noi,nof);
for iter_no = 1:noi
    tic;
    for i=nof-1:-1:1
        [Image_reg(:,:,i),shifts(:,iter_no,i)] = Reg_GS_rigid(Image(:,:,i),mean(Image_reg(:,:,i+1:end),3),step,10,mask);
    end
    Image = Image_reg;
    step = step/2;
    fprintf(['Iteration = ' num2str(iter_no) '...']);
    toc;
end