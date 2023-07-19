function temp = low_rank(Image,bloc_x,bloc_y,sx,nof,tau)


temp = reshape(Image,[bloc_x,sx/bloc_x,bloc_y,sx/bloc_y,nof]);
temp = permute(temp,[1 3 2 4 5]);
temp = reshape(temp,[bloc_x*bloc_y,sr^2/bloc_x/bloc_y,nof]);

if(bloc_x*bloc_y<nof)
    temp = permute(temp,[3 1 2]);
else
    temp = permute(temp,[1 3 2]);
end
temp = gather(temp); % CPU is fater here. When using parfor, ~10 times faster than old code.

parfor i=1:sx^2/bloc_x/bloc_y
    [U,S,V] = givefastSVD(temp(:,:,i));
    S = diag(S)-tau;
    S(S<0) = 0;
    temp(:,:,i) = U*diag(S)*V';
end

if(bloc_x*bloc_y<nof)
    temp = permute(temp,[2 3 1]);
else
    temp = permute(temp,[3 2 1]);
end

temp = gpuArray(reshape(temp,[bloc_x,bloc_y,sx/bloc_x,sx/bloc_y,nof]));
temp = permute(temp,[1 3 2 4 5]);
temp = reshape(temp,[sr,sr,sz]);

%Image = Image + lambda_new * (lr_term_update-Image);   
%clear lr_term_update
end

function [u,s,v] = givefastSVD(X)

[v,s] = eig(X'*X);
s = sqrt(abs(s));

s = diag(s);
index = s >1e-6*max(s);
v = v(:,index);
s = diag(s(index));
u = X*v;

u = bsxfun(@rdivide,u,sqrt(sum(abs(u).^2,1)));
end
