function T1_map = T1_fit_MOLLI(im,InversionTime)

siz = size(im);
im = reshape(im,siz(1)*siz(2),siz(3));


im = repmat(im,[1,1,4]);
im(:,1,1:4) = -im(:,1,1:4);
im(:,2,2:4) = -im(:,2,2:4);
im(:,3,3:4) = -im(:,3,3:4);
im(:,4,4) = -im(:,4,4);

f = @(A,B,T1_star,x) A-B.*exp(-x./T1_star);


for iset=1:4
    parfor i=1:siz(1)*siz(2)
        i
        [A_temp(i,iset),B_temp(i,iset),T1_temp(i,iset)] = T1_fitting(InversionTime',im(i,:,iset)',f);
    end
end

mse = A_temp- B_temp.*exp(-reshape(InversionTime,[1,1,1,siz(3)])./T1_temp);
mse = squeeze(sum((permute(mse,[1,4,2,3]) - im).^2,2));
[~,idx] = min(mse,[],2);

for i=1:siz(1)*siz(2)
    A_final(i) = A_temp(i,idx(i));
    B_final(i) = B_temp(i,idx(i));
    T1_final(i) = T1_temp(i,idx(i));
end
% 
% A_final = zeros(sx,sy);
% B_final = zeros(sx,sy);
% T1_star_final = zeros(sx,sy);
% 
% for i=1:sx
%     for j=1:sy
%         A_final(i,j) = A(i,j,idx(i,j));
%         B_final(i,j) = B(i,j,idx(i,j));
%         T1_star_final(i,j) = T1_star(i,j,idx(i,j));
%     end
% end

T1_map = T1_final.*(B_final./A_final - 1);
T1_map = reshape(T1_map,siz(1),siz(2));
T1_map = abs(T1_map);
T1_map(T1_map>4) = 4;

end


function [A,B,T1_star] = T1_fitting(t,intensity,f)

warning off
t1_fit = fit(t,intensity,f,'StartPoint',[0,0,50],'Upper',[1e5,1e5,1e7],'Lower',[0,0,0]);

A = t1_fit.A;
B = t1_fit.B;
T1_star = t1_fit.T1_star;
end
