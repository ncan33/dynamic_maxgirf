function out = SSIM(im1,im2)

im1 = double(im1);
im2 = double(im2);

for i=1:size(im1,3)
    for j=1:size(im1,4)
        ss(1,i,j) = ssim(im1(:,:,i,j),im2(:,:,i,j));
    end
end

out = mean(mean(ss,2),3);
