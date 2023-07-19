function kSpace = exact_fft2(image, kx, ky)

[sx, sy, nof, nc] = size(image);
skx = size(kx,1);
[x, y] = meshgrid(1:sx,1:sy);
nor = size(kx,2);
kcenter = sx/2;
x = x - kcenter;
y = y - kcenter;
if isa(image,'gpuArray')
    kSpace = gpuArray(single(zeros(sx,nor,nof,nc)));
else
    kSpace = single(zeros(sx,nor,nof,nc));
end
%{
for Ncoil = 1:nc
    Ncoil
    for Nframe = 1:nof
        Nframe
        for Nx = 1:sx
        for Nrays = 1:nor
            kx_temp = kx(Nx,Nrays,Nframe);
            ky_temp = ky(Nx,Nrays,Nframe);
            temp = image(:,:,Nframe,Ncoil) .* exp(-2*pi*1i*(y*kx_temp/sx+x*ky_temp/sx));
            kSpace(Nx,Nrays,Nframe,Ncoil) = sum(sum(temp));
        end
        end
    end
end
%}

%for Ncoil = 1:nc

%    for Nframe = 1:nof
        tic
        for Nx = 1:skx
            %Nx
        for Nrays = 1:nor
            kx_temp = kx(Nx,Nrays,:);
            ky_temp = ky(Nx,Nrays,:);
            temp = bsxfun(@times,image,exp(-2*pi*1i*(bsxfun(@times,y,kx_temp)/skx+bsxfun(@times,x,ky_temp)/skx)));
            kSpace(Nx,Nrays,:,:) = sum(sum(temp,1),2);
        end
        end
        toc
%    end
%end
