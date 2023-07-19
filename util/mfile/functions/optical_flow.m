function dd = optical_flow(im1, im2, para)

[sx, sy] = size(im1);

c1 = ones(sx, sy, class(im1));
c2 = ones(sx, sy, class(im1));

for ilevel = 1:para.levels + 1
    im1_temp = im1;
    im2_temp = im2;
    c1_temp = c1;
    c2_temp = c2;
    
    if ilevel ~= para.levels + 1
        for i = 1:para.levels + 1 - ilevel
            im1_temp = impyramid(im1_temp, 'reduce');
            im2_temp = impyramid(im2_temp, 'reduce');
            c1_temp = impyramid(c1_temp, 'reduce');
            c2_temp = impyramid(c2_temp, 'reduce');
        end
    end
    
    if ilevel == 1
        dx = zeros(size(im1_temp), class(im1_temp));
        dy = zeros(size(im1_temp), class(im1_temp));
        dd = cat(3, dx, dy);
    end
    
    if size(dd, 1) ~= size(im1_temp, 1) | size(dd, 2) ~= size(im1_temp, 2)
        dd = imresize(dd, size(im1_temp));
    end

    dd = cal_flow_iteration(im1_temp, im2_temp, c1_temp, c2_temp, dd * 2, para);

    if ilevel ~= para.levels + 1
        dd = impyramid(dd, 'expand');
    end

end

end



