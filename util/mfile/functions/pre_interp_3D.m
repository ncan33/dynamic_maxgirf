function [Data,para] = pre_interp_3D(kSpace_radial,kx,ky,para)
disp('Pre-interpolate into Cartesian space...');t1=tic;
%keyboard
%interp_dir = para.Recon.cart_kSpace_dir;
%ifInterp   = para.Recon.ifInterp;

%if ifInterp == 1
%    load(interp_dir);
%    para.CPUtime.prepare_kSpace = toc;toc;fprintf('\n');
%    return
%end

switch para.Recon.interp_method
    case 'GROG'
        Data.kSpace = GROG.GROG_3D( kSpace_radial, kx, ky, 0, 1);
    case 'Toeplitz'
        [Data.kSpace,Data.mask,Data.Apodizer] = Toeplitz_3D(kSpace_radial,kx,ky,para);
    case 'NUFFT'
        Data = NUFFT.ThreeD_init(kSpace_radial,kx,ky,para);
    case 'grid3'
        [sx,nor,sz,nof,no_comp] = size(kSpace_radial);
        data = double(reshape(kSpace_radial,[sx*nor,sz*nof*no_comp]));
        kx = double(kx);
        ky = double(ky);
        
        x = reshape(kx,[sx*nor,sz*nof]);
        y = reshape(ky,[sx*nor,sz*nof]);
        
        x = repmat(x,[1 no_comp]);
        y = repmat(y,[1 no_comp]);

        Xr = round(x);
        Yr = round(y);
        
        kSpace_cart = single(zeros((sx+1)*(sx+1),sz*nof*no_comp));
        kSpace_r = single(zeros(sx*nor,sz*nof*no_comp));
        
        for i=1:sz*nof*no_comp
            warning off
            index = data(:,i) ~= 0;
            kSpace_r(index,i) = griddata(x(index,i),y(index,i),data(index,i),Xr(index,i),Yr(index,i));
        end
        
        kSpace_r(isnan(kSpace_r)) = 0;
        
        indx = sub2ind([sx+1,sx+1,sz*nof*no_comp],Xr+sx/2+1,Yr+sx/2+1);
        
        for i=1:sz*nof*no_comp
            kSpace_cart(indx(:,i),i) = kSpace_r(:,i);
        end
        
        kSpace_cart = reshape(kSpace_cart,[sx+1,sx+1,sz,nof,no_comp]);
        kSpace_cart(1,:,:,:,:,:) = [];
        kSpace_cart(:,1,:,:,:,:) = [];
        Data.kSpace = kSpace_cart;
end
%save(interp_dir,'kSpace_cart','-v7.3')
para.CPUtime.prepare_kSpace = toc(t1);toc(t1);fprintf('\n');
end