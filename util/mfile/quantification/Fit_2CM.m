function [param,model_curves] = Fit_2CM(AIF,Tissue)
% param = Fit_2CM(AIF,Tissue)
K_trans = 0.8;
K_ep= 2;
vb = 0.01;
dt = 0.05;
t = (0:1:size(Tissue,2)-1)/60;

in = {double(AIF),double(t)};

para = [K_trans,K_ep,vb,dt];
para_low = [0,0,0,-0.1];
para_up = [10,30,1,0.1];

options = optimoptions('lsqcurvefit','Display','off');
for i=1:size(Tissue,1)
    param(:,i) = lsqcurvefit(@Two_component_model,para,in,double(Tissue(i,:)),para_low,para_up,options);
    model_curves(i,:) = Two_component_model(param(:,i),in);
end

figure(1)
clf
hold on
col = colors;
col = [col;col/2];
for i=1:size(Tissue,1)
    plot(Two_component_model(param(:,i),in),'Color',col(i,:))
    plot(Tissue(i,:),'Color',col(i,:),'Marker','.','LineStyle','none')
end

end

function c = Two_component_model(x,y)

AIF = y{1};
t = y{2};

K_trans = x(1);
K_ep = x(2);
vb = abs(x(3));
dt = x(4);

AIF = interp1((0:length(AIF)-1)/60,AIF,t-dt);
AIF(isnan(AIF)) = 0;


%AIF = [zeros(1,dt),AIF];
N = length(t);
c = K_trans*conv(AIF,exp(-K_ep*t));
c = c(1:N);
c = c/60 + vb*AIF(1:N);
c = c;

end