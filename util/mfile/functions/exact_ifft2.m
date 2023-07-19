function image = exact_ifft2(kSpace,kx,ky)

[sx,nor,nof,nc] = size(kSpace);
image = zeros(sx,sx,nof,nc);
kcenter = sx/2 - 1;

for Ncoil = 1:nc
    for Nframe = 1:nof
        kx_temp = kx(:,:,Nframe);
        ky_temp = ky(:,:,Nframe);
        %k1 = bsxfun(@times,kSpace_87(:,:,1,Ncoil),W);
        for Nx = 1:sx
            for Ny = 1:sx
                x = Nx-kcenter;
                y = Ny-kcenter;
                temp = kSpace(:,:,Nframe,Ncoil) .* exp(2*pi*1i*(x*kx_temp/sx+y*ky_temp/sx));
                image(Nx,Ny,Ncoil) = sum(sum(temp));
            end
        end
    end
end
% imagesc(abs(image))