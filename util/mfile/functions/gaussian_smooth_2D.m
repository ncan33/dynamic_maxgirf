function out = gaussian_smooth_2D(Image,sigma)

Image = abs(Image);
Image = medfilt2(Image,[5,5]);
[sx,sy] = size(Image);
out = Image;
localmax = imregionalmax(Image);
idx = find(localmax);
idx(:,2) = Image(idx);
idx = sortrows(idx,2,'descend');
%{
for i=1:200
    idx_temp = idx(i,1);
    [x_temp,y_temp] = ind2sub([sx,sy],idx_temp);
    [Y,X] = meshgrid(1:sx,1:sy);
    Y = Y - y_temp;
    X = X - x_temp;
    G = idx(i,2)*exp(-(X.^2+Y.^2)/sigma^2);
    mask = (out-G)<0;
    out(mask) = G(mask);
end
%}
out = out/max(out(:));
out(out<0.1) = 0.1;
for n=1:3
for i=1:sx
    temp = out(i,:);
    temp_smooth = smooth(temp,0.1);
    temp = max(temp,temp_smooth.');
    out(i,:) = temp;
end
for i=1:sy
    temp = out(:,i);
    temp_smooth = smooth(temp,0.1);
    temp = max(temp,temp_smooth);
    out(:,i) = temp;
end
out = imgaussfilt(out,1.5);
end
