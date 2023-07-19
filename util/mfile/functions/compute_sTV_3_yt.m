function stv = compute_sTV_3_yt(img,weight,beta_sqrd)
siz = size(img);

diffx = diff(img,1,1);
diffy = diff(img,1,2);

d1 = cat(1,zeros([1 siz(2:end)]),diffx);
d2 = cat(2,zeros([siz(1) 1 siz(3:end)]),diffy);

d3 = cat(1,diffx,zeros([1 siz(2:end)]));
d4 = cat(2,diffy,zeros([siz(1) 1 siz(3:end)]));

d5 = d2; d5(1,:,:,:,:,:) = []; d5 = cat(1,d5,zeros([1,siz(2:end)]));
d6 = d1; d6(:,1,:,:,:,:) = []; d6 = cat(2,d6,zeros([siz(1),1,siz(3:end)]));


t1 = (d1 + d2)./sqrt(abs(d1).^2 + abs(d2).^2 + beta_sqrd);

t2 = - d3./sqrt(abs(d3).^2 + abs(d5).^2 + beta_sqrd);

t3 = - d4./sqrt(abs(d4).^2 + abs(d6).^2 + beta_sqrd);

stv = weight*(t1+t2+t3);
