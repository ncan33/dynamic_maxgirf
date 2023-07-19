function [I2,shifts,E] = Reg_GS_rigid(I1,I0,step,noi,mask)
[sx,sy] = size(I0);
Maxiter = noi;
a = 30;
s = 3;

kx = cos(2*pi*(0:sx/(sx-1):sx));
ky = cos(2*pi*(0:sy/(sy-1):sy));
W = 2*(kx+ky.'-2);
W = single((1-a*W).^-s);
 
%[y,x] = meshgrid(1:sx,1:sy);

I2_fft = fftshift2(fft2(fftshift2(I1)));

PhaseX = reshape(-2i*pi*((1:sx)-floor(sx/2)-1)/sx,sx,1);
PhaseY = reshape(-2i*pi*((1:sy)-floor(sy/2)-1)/sy,1,sy);
Phase = exp(PhaseX*0 + PhaseY*0);
shifts = zeros(2,1);
for iter=1:Maxiter
    
    I2 = fftshift2(ifft2(fftshift2(I2_fft.*Phase)));
    I2 = abs(I2);
    
    ddx = 0.5*(I2(3:end,:) - I2(1:end-2,:));
    ddy = 0.5*(I2(:,3:end) - I2(:,1:end-2));
    ddx = cat(1,I2(2,:) - I2(1,:),ddx,I2(end,:) - I2(end-1,:));
    ddy = cat(2,I2(:,2) - I2(:,1),ddy,I2(:,end) - I2(:,end-1));
    
    dI = I2-I0;
    %dI = dI - medfilt2(dI,[3,3]);

    %J = jacobian(y,x);

    %dx = W.*fft2(ddx.*dI.*J);
    dx = W.*fft2(ddx.*dI);
    dx = -real(ifft2(dx));

    %dy = W.*fft2(ddy.*dI.*J);
    dy = W.*fft2(ddy.*dI);
    dy = -real(ifft2(dy));
    
    dx = dx.*mask;
    dy = dy.*mask;

    d = sqrt(dx.^2+dy.^2);

    md = max(d(:));
 
    dx = dx/md;
    dy = dy/md;
    
    dx = mean(mean(-dx*step));
    dy = mean(mean(-dy*step));
    
    E(iter) = sum(d(:));
    if iter>1 && E(iter)>E(iter-1)
        step = step/2;
        %return
    end
%     
%     figure(2)
%     plot(E)
%     drawnow
    
    
    Phase = exp(PhaseX*dx + PhaseY*dy);

    shifts = shifts + [dx;dy];

    I2_fft = fftshift2(fft2(fftshift2(I2)));
    
%     x = x + step*dx;
%     y = y + step*dy;
%     
%     x(x<1) = 1;
%     x(x>sx) = sx;
%     y(y<1) = 1;
%     y(y>sy) = sy;

%     figure(1)
%     subplot(2,1,1)
%     imagesc([I2,I0,dI])
%     colormap gray
%     axis image
%     brighten(0.4)
%     drawnow

end
end
function J = jacobian(x,y)
    %[Y,X] = meshgrid(1:sx,1:sy);
    [fxx,fxy] = gradient(x);
    [fyx,fyy] = gradient(y);
    J = fxx .* fyy - fxy .* fyx;
end

function output = fftshift2(input)

output = fftshift(input,1);
output = fftshift(output,2);

end

