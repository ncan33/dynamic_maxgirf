function update = low_rank_multiple_coil(Image, bloc_x, bloc_y, tau)

[sx, sy, nc, nof] = size(Image);
n_block_x = sx/bloc_x;
n_block_y = sy/bloc_y;

update = reshape(Image,[bloc_x, n_block_x, bloc_y, n_block_y, nc, nof]);
update = permute(update,[1, 3, 5, 2, 4, 6]);
update = reshape(update,[bloc_x * bloc_y * nc, n_block_x * n_block_y, nof]);
update = permute(update,[3, 1, 2]);
% update = gather(update); 
f = @(x) svd(x, 0);
%     [U,S,V] = givefastSVD(update(:,:,i));
%     [U,S,V] = svd(update(:,:,i),0);

% CPU is fater here. When using parfor, ~10 times faster than old code.
parfor i = 1:n_block_x * n_block_y

    [U,S,V] = f(update(:,:,i));

    S = S - tau;
    S(S<0) = 0;
%     fprintf([num2str(max(S(:))),'\n'])
    update(:, :, i) = U * S * V';
end

update = permute(update, [2, 3, 1]);
% update = gpuArray(reshape(update,[bloc_x,bloc_y,sx/bloc_x,sy/bloc_y,nof]));
update = reshape(update, [bloc_x, bloc_y, nc, n_block_x, n_block_y, nof]);
update = permute(update, [1, 4, 2, 5, 3, 6]);
update = reshape(update, [sx, sy, nc, nof]);
end