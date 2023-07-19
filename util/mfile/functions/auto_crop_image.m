function Image = auto_crop_image(Image,siz)

[sx,sy,sz,nof] = size(Image);

std_map = std(Image,1,4);
std_smoothed = imgaussfilt(std_map,3);
std_smoothed = reshape(std_smoothed,sx*sy,sz);
std_smoothed = sort(std_smoothed,'descend');
threshold = std_smoothed(round(sx*sy*0.01),:);
mask = false([sx,sy,sz]);
for i=1:sz
    mask(:,:,i) = std_map(:,:,i) >= threshold(i);
end



box = ones(siz,'single');
Image = abs(Image);
std_map = std(Image,1,4);
d_map = squeeze(max(Image,[],4) - min(Image,[],4));
nof = size(d_map,3);
d_map_conv = zeros(size(Image,1),size(Image,2),size(Image,3));
for i=1:nof
    d_map_conv(:,:,i) = conv2(d_map(:,:,i),box,'same');
end
d_map_conv = sum(d_map_conv,3);
localmax = imregionalmax(d_map_conv);
[y,x] = find(localmax);
distance = sos([x,y] - [size(d_map_conv,1),size(d_map_conv,2)]/2);
[~,idx] = min(distance);
x = x(idx);
y = y(idx);

% x = round(x-siz(1)/2);
% y = round(y-siz(2)/2);

x = x + [1:siz(1)] - round(siz(1)/2);
y = y + [1:siz(2)] - round(siz(2)/2);

Image_crop = Image(x,y,:,:);
keyboard


