function res = im2row3d(im, winSize)
%res = im2row(im, winSize)
[sx,sy,sz,ncoil] = size(im);

res = zeros((sx-winSize(1)+1)*(sy-winSize(2)+1)*(sz-winSize(3)+1), prod(winSize), ncoil, 'like', im);
count = 0;
for y = 1:winSize(2)
    for x = 1:winSize(1)
        for z = 1:winSize(3)
            count = count + 1;
            res(:,count,:) = reshape(im(x:sx-winSize(1)+x, y:sy-winSize(2)+y, z:sz-winSize(3)+z, :),...
                (sx-winSize(1)+1)*(sy-winSize(2)+1)*(sz-winSize(3)+1), ncoil);
        end
    end
end
