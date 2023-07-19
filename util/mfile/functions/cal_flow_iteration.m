function d = cal_flow_iteration(im1, im2, c1, c2, displacement, para)

sigma_flow = para.sigma_flow;
iterations = para.iterations;
sigma_poly = para.sigma_poly;

[sx, sy] = size(im1);

dx = displacement(:, :, 1);
dy = displacement(:, :, 2);

[A1, B1, C1] = poly_exp(im1, c1, sigma_poly);
[A2, B2, C2] = poly_exp(im2, c2, sigma_poly);

[y, x] = meshgrid(1:sx, 1:sy);

n_flow = round(4 * sigma_flow + 1);
xw = -n_flow : n_flow;
w = exp(-xw.^2 / (2 * sigma_flow.^2));

S = eye(2);
S_T = S';

for iiter = 1:iterations
   dx_ = round(dx);
   dy_ = round(dy);
   
   x_ = x + dx_;
   y_ = y + dy_;
   
   x_2 = max(min(x_, sx), 1);
   y_2 = max(min(y_, sy), 1);
   
   off_x = logical(x_ - x_2);
   off_y = logical(y_ - y_2);
   off_ = off_x | off_y;
   
   x_ = x_2;
   y_ = y_2;
   
   for i = 1:sx
       for j = 1:sy
           c_(i, j) = c1(x_(i, j), y_(i, j));
       end
   end
   
   c_(off_) = 0;
   
   for i = 1:sx
       for j = 1:sy
           A(i, j, :, :) = A1(i, j, :, :) + A2(x_(i, j), y_(i, j), :, :);
       end
   end
   A = A/2;
   
   A = A .* c_;
   
   for i = 1:sx
       for j = 1:sy
           delB(i, j, :) = -0.5 * (B2(x_(i, j), y_(i, j), :) - B1(i, j, :));
       end
   end
   delB = delB + sum(A .* cat(4, dx_, dy_), 4);
   delB = delB .* c_;
   
   A_T = permute(A, [1, 2, 4, 3]);
   ATA = squeeze(sum(permute(S_T, [3, 4, 1, 2]) .* permute(A_T, [1, 2, 5, 3, 4]), 4));
   ATA = squeeze(sum(ATA .* permute(A, [1, 2, 5, 3, 4]), 4));
   ATA = squeeze(sum(ATA .* permute(S, [3, 4, 5, 1, 2]), 4));
   
   ATb = squeeze(sum(permute(S_T, [3, 4, 1, 2]) .* permute(A_T, [1, 2, 5, 3, 4]), 4));
   ATb = squeeze(sum(ATb .* permute(delB, [1, 2, 4, 3]), 4));
   
   G_avg = squeeze(mean(mean(ATA)));
   h_avg = squeeze(mean(mean(ATb)));
   p_avg = linsolve(G_avg, h_avg);
   d_avg = S * p_avg;
   
   if iiter == 1
       mu = 0.5 * trace(G_avg);
   end
   
   A_TA = squeeze(sum(A_T .* permute(A, [1, 2, 5, 3, 4]), 4));
   A_TdelB = squeeze(sum(A_T .* permute(delB, [1, 2, 4, 3]), 4));
   for i = 1:size(A_TA, 3)
       for j = 1:size(A_TA, 4)
           G(:, :, i, j) = conv2(A_TA(:, :, i, j), flip(w)', 'same');
           G(:, :, i, j) = conv2(G(:, :, i, j), flip(w), 'same');
       end
       h(:, :, i) = conv2(A_TdelB(:, :, i), flip(w)', 'same');
       h(:, :, i) = conv2(h(:, :, i), flip(w), 'same');
   end
   
   G = G + permute(mu * eye(2), [3, 4, 1, 2]);
   h = h + permute(mu * d_avg, [3, 4, 1, 2]);
   
   for i = 1:sx
       for j = 1:sy
           d(i, j, :) = linsolve(squeeze(G(i, j,:, :)), squeeze(h(i, j, :)));
       end
   end
   
   dx = d(:, :, 1);
   dy = d(:, :, 2);
   
end


end