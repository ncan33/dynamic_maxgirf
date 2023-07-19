function tTV_update = compute_tTV_pixel_interp(Image,weight,beta_square,Motion)

nof  = size(Image,3);
ns   = size(Image,4);
nsms = size(Image,5);

for j=1:ns
    for k=1:nsms
        for i=2:nof-1
            I0 = Image(:,:,i,j,k);
            If = Image(:,:,i+1,j,k);
            Ib = Image(:,:,i-1,j,k);
            
            If = interp2(If,Motion.xb(:,:,i,j,k),Motion.yb(:,:,i,j,k),'cubic');
            Ib = interp2(Ib,Motion.xf(:,:,i,j,k),Motion.yf(:,:,i,j,k),'cubic');
            
            temp_f = If - I0;
            temp_b = I0 - Ib;
            
            temp_f = temp_f./sqrt(abs(temp_f).^2+beta_square);
            temp_b = temp_b./sqrt(abs(temp_b).^2+beta_square);
            
            tTV_update(:,:,i,j,k) = temp_f-temp_b;
        end
        
        tTV_update(:,:,1,j,k) = interp2(Image(:,:,2,j,k),Motion.xb(:,:,1,j,k),Motion.yb(:,:,1,j,k),'cubic') - Image(:,:,1,j,k);
        tTV_update(:,:,nof,j,k) = interp2(Image(:,:,end-1,j,k),Motion.xf(:,:,end,j,k),Motion.yf(:,:,end,j,k),'cubic') - Image(:,:,end,j,k);
    end
end
%tTV_update(:,:,end+1,:,:) = 0;
%tTV_update(:,:,1,:,:) = temp_f(:,:,1,:,:);
tTV_update = weight*tTV_update;
