function im_fill = spline_2D_yt(im)

%mask = im>mean(im(:));
im_smooth = imfilter(im,fspecial('gaussian',size(im),3));
mask = im_smooth>mean(im(:));
im_smooth = imfilter(im.*mask,fspecial('gaussian',size(im),1));
im_regionalmax = imregionalmax(im_smooth);
[x,y] = find(im_regionalmax);
%figure,imagesc(im_smooth);
%hold on
%plot(y,x,'ro')
xy = [y.';x.'];
values = im(im_regionalmax);
p = tpaps(xy,values.',0.0002);

[X,Y] = meshgrid(1:size(im,1),1:size(im,2));
XY = [X(:).';Y(:).'];
im_spline = fnval(p,XY);
im_spline = reshape(im_spline,size(im));

mean_im = mean(im(:));
im_fill = regionfill(im_spline,im_spline<mean_im);

d = im_fill - im_spline;
if sum(d(:)) == 0
    return
end
max_d = max(d(:));
im_fill = im_spline + d/max_d*(max_d-mean_im/2);
