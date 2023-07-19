function [Im_reg, x, y] = registration_2D_to_one_frame(Image,frame)

%Image_reg = imgaussfilt(abs(Image),1);
Im_reg = zeros(size(Image),'like',Image);
x = zeros(size(Image),'like',Image);
y = zeros(size(Image),'like',Image);
Im_reg(:,:,frame) = Image(:,:,frame);
nof = size(Image,3);

for i = frame+1:nof
    [Im_reg(:,:,i), x(:,:,i), y(:,:,i)] = Reg_GS_accurate(Image(:,:,i),Im_reg(:,:,i-1),0.1,100);
    %Im_reg(:,:,i) = interp2(Image(:,:,i),y,x);
end

for i = frame-1:-1:1
    [Im_reg(:,:,i), x(:,:,i), y(:,:,i)] = Reg_GS_accurate(Image(:,:,i),Im_reg(:,:,i+1),0.1,100);
    %Im_reg(:,:,i) = interp2(Image(:,:,i),y,x);
end