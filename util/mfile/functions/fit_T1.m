d_all = 0;
for i=1:20
    temp = Image_process(:,:,:,i);
    d = sum((temp-M0).^2,3);
    [d,order] = min(d,[],4);
    T1_temp(:,:,i) = T1(order);

end
