function idx = find_line_spiral(kx,ky,r)

r = 84;
dtheta = 1/r;
n_line = round(pi/dtheta);

idx = zeros(n_line,r*2+1);

figure
plot(kx,ky,'.')
axis image
hold on

for i=1:n_line %:r-1
    line = -r:r;
    theta = dtheta*i;
    x = line*cos(theta);
    y = line*sin(theta);
    plot(x,y);
    for j=1:length(line)
        
        x_temp = x(j);
        y_temp = y(j);
        d = sqrt((kx-x_temp).^2 + (ky-y_temp).^2);
        [d_min,idx_temp] = min(d(:));
%         if d_min<1
            idx(i,j) = idx_temp;
%         else
%             d_min
%         end
    end
end
