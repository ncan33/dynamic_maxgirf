function show_SMS_yt_pause(img, brightness)

if nargin < 2
    brightness = 0.2;
end
if ~isreal(img)
    img = abs(img);
end
[sx,sy,nof,nSMS,ns] = size(img);
img = reshape(img,[sx sy nof nSMS*ns]);

img = permute(img,[1 2 4 3]);

img = reshape(img,[sx,sy*nSMS*ns,nof]);

figure
for i=1:size(img,3)
    im_temp = img(:,:,i);
    
    imagesc(im_temp)
    colormap gray
    brighten(brightness)
    axis equal
    title(i)

    pause
end
end