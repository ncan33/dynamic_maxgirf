function X = fft3GPU(X)
keyboard
siz = size(X);
X = reshape(X,[siz(1:3),prod(siz(4:end))]);
for i=1:prod(siz(4:end))
    X(:,:,:,i) = fftn(X(:,:,:,i));
end
X = reshape(X,siz);
end