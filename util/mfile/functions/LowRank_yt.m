function Image = LowRank_yt(Image)
siz = size(Image);
if length(siz) < 4
    siz(4) = 1;
end
Image = reshape(Image,[prod(siz(1:2)),siz(3),prod(siz(4:end))]);


%N = siz(3) - round(siz(3)*0.1);
N = siz(3);

for i=1:prod(siz(4:end))
    Image(:,:,i) = low_rank_update(Image(:,:,i),N);
end

Image = reshape(Image,siz);
end

function Image = low_rank_update(Image,N)
[U,S,V] = svd(Image,0);
S = diag(S);
S = S - S(1)*0.02;
S(S<0) = 0;
S = diag(S);
Image = U*S*V';
end