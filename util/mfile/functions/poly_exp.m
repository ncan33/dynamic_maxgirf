function [A, B, C] = poly_exp(im, c, sigma_poly)
% im = ones(13);
% c = c1;
% sigma_poly = 0.5;

n = round(4 * sigma_poly + 1);
x = -n : n;
a = exp(-x.^2 / (2 * sigma_poly.^2));

bx = [ones(size(a)); x; ones(size(a)); x.^2; ones(size(a)); x]';
by = [ones(size(a)); ones(size(a)); x; ones(size(a)); x.^2; x]';

cf = c .* im;

G = zeros([size(im), size(bx, 2), size(bx, 2)]);
v = zeros([size(im), size(bx, 2)]);

ab = a' .* bx;
abb = ab .* permute(bx, [1, 3, 2]);

for i = 1:size(abb, 2)
    for j = 1:size(abb, 3)
        G(:, :, i, j) = conv2(c, flip(abb(:, i, j)), 'same');
    end
    v(:, :, i) = conv2(cf, flip(ab(:, i)), 'same');
end

ab = a' .* by;
abb = ab .* permute(by, [1, 3, 2]);

for i = 1:size(abb, 2)
    for j = 1:size(abb, 3)
        G(:, :, i, j) = conv2(G(:, :, i, j), flip(abb(:, i, j))', 'same');
    end
    v(:, :, i) = conv2(v(:, :, i), flip(ab(:, i))', 'same');
end

for i = 1:size(G, 1)
    for j = 1:size(G, 2)
        r(i, j, :, :)  = linsolve(squeeze(G(i, j, :, :)), squeeze(v(i, j, :)));
    end
end

A = r(:, :, 4);
A(:, :, 1, 2) = r(:, :, 6) / 2;
A(:, :, 2, 1) = A(:, :, 1, 2);
A(:, :, 2, 2) = r(:, :, 5);

B = r(:, :, 2);
B(:, :, 2) = r(:, :, 3);

C = r(:, :, 1);
end
