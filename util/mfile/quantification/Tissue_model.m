function c_all = Tissue_model(AIF,x)

t = x{2};
x = x{1};
AIF = Quant.AIF_model(AIF,t);

N = size(x,2);
K_trans = x(1,:);
K_ep = x(2,:);
vb = abs(x(3,:));
dt = x(4,:);

t = t/60;

for i=1:N
    AIF_temp = interp1((0:length(AIF)-1)/60,AIF,t-dt(i));
    AIF_temp(isnan(AIF_temp)) = 0;
    
    n = length(t);
    c = K_trans(i)*conv(AIF_temp,exp(-K_ep(i)*t));
    c = c(1:n);
    c = c/60 + vb(i)*AIF_temp(1:n);
    c_all(:,i) = c;
end

