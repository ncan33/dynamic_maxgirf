function [Image_reg] = deformable_reg_GS(Image)

% if ~exist('mask')
%     mask = get_mask(Image(:,:,end));
% end
Image_reg = Image;
nof = size(Image,3);
noi = 4;
step = 1;

%shifts = zeros(2,noi,nof);
% for iter_no = 1:noi
%     tic;
%     for i=nof-1:-1:1
%         [Image_reg(:,:,i)] = Reg_GS_tv(Image(:,:,i),mean(Image_reg(:,:,i+1:end),3),step,10);
%     end
%     Image = Image_reg;
%     step = step/2;
%     fprintf(['Iteration = ' num2str(iter_no) '...']);
%     toc;
% end


for iter_no = 1:noi
    tic;
    for i=nof-2:-1:3
        [Image_reg(:,:,i)] = Reg_GS_accurate(Image(:,:,i),mean(Image(:,:,i-2:i+2),3),step,10);
    end
    Image = Image_reg;
    step = step/2;
    fprintf(['Iteration = ' num2str(iter_no) '...']);
    toc;
end