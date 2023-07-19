function [order,order_back] = get_order(Image)

[sx,sy,nof] = size(Image);

[y,x] = meshgrid(1:sx,1:sy);
x = repmat(x,[1,1,nof]);
y = repmat(y,[1,1,nof]);

[~,order] = sort(abs(Image),3);
[~,order_back] = sort(order,3);

order = sub2ind([sx,sy,nof],x,y,order);
order_back = sub2ind([sx,sy,nof],x,y,order_back);