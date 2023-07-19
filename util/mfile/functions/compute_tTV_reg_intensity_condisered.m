function tTV_update = compute_tTV_reg_intensity_condisered(image,weight,beta_square)
%tTV_update = compute_tTV_yt(image,weight,beta_square)
if weight~=0
    temp_a = diff(image,1,3);
    for i=1:size(temp_a,3)
        for j=1:size(temp_a,5)
            temp_a(:,:,i,1,j) = temp_a(:,:,i,1,j) - medfilt2(real(temp_a(:,:,i,1,j)),[9,9]) - 1i*medfilt2(imag(temp_a(:,:,i,1,j)),[9,9]);
        end
    end
    temp_b = temp_a./(sqrt(beta_square+(abs(temp_a).^2)));
    temp_c = diff(temp_b,1,3);
    tTV_update = weight * cat(3,temp_b(:,:,1,:,:,:),temp_c,-temp_b(:,:,end,:,:,:));
else
    tTV_update = 0;
end
return

[s1,s2,s3,s4,s5] = size(image);
image = image(:,:,:,:);

im_real = real(image);
im_imag = imag(image);

nof = s3;
ns  = s4*s5;

for i = 1:ns
    for j = 1:nof-1
        Ib_real = im_real(:,:,j,i);
        If_real = im_real(:,:,j+1,i);
        
        Ib_imag = im_imag(:,:,j,i);
        If_imag = im_imag(:,:,j+1,i);
        
        Ib_real_all(:,:,j,i) = Reg_GS_tv(If_real,Ib_real,0.5,20);
        If_real_all(:,:,j,i) = Reg_GS_tv(Ib_real,If_real,0.5,20);
        
        Ib_imag_all(:,:,j,i) = Reg_GS_tv(If_imag,Ib_imag,0.5,20);
        If_imag_all(:,:,j,i) = Reg_GS_tv(Ib_imag,If_imag,0.5,20);
        
    end
end

d_real_lmr = Ib_real_all - im_real(:,:,1:end-1,:);
d_real_rml = im_real(:,:,2:end,:) - If_real_all;
d_imag_lmr = Ib_imag_all - im_imag(:,:,1:end-1,:);
d_imag_rml = im_imag(:,:,2:end,:) - If_imag_all;

for i = 1:ns
    for j = 1:nof-1
        d_real_lmr(:,:,j,i) = d_real_lmr(:,:,j,i) - medfilt2(d_real_lmr(:,:,j,i),[3,3]);
        d_real_rml(:,:,j,i) = d_real_rml(:,:,j,i) - medfilt2(d_real_rml(:,:,j,i),[3,3]);
        d_imag_lmr(:,:,j,i) = d_imag_lmr(:,:,j,i) - medfilt2(d_imag_lmr(:,:,j,i),[3,3]);
        d_imag_rml(:,:,j,i) = d_imag_rml(:,:,j,i) - medfilt2(d_imag_rml(:,:,j,i),[3,3]);
    end
end

d_real_lmr = d_real_lmr./sqrt(beta_square+d_real_lmr.^2);
d_real_rml = d_real_rml./sqrt(beta_square+d_real_rml.^2);
d_imag_lmr = d_imag_lmr./sqrt(beta_square+d_imag_lmr.^2);
d_imag_rml = d_imag_rml./sqrt(beta_square+d_imag_rml.^2);

update_real = d_real_lmr(:,:,2:nof-1,:) - d_real_rml(:,:,1:nof-2,:);
update_imag = d_imag_lmr(:,:,2:nof-1,:) - d_imag_rml(:,:,1:nof-2,:);

update_real = cat(3,d_real_lmr(:,:,1,:),update_real,-d_real_rml(:,:,end,:));
update_imag = cat(3,d_imag_lmr(:,:,1,:),update_imag,-d_imag_rml(:,:,end,:));

tTV_update = update_real + 1i*update_imag;
tTV_update = weight * reshape(tTV_update,[s1,s2,s3,s4,s5]);
return

    
if weight~=0
    temp_a = diff(image,1,3);
    temp_b = temp_a./(sqrt(beta_square+(abs(temp_a).^2)));
    temp_c = diff(temp_b,1,3);
    tTV_update = weight * cat(3,temp_b(:,:,1,:,:,:),temp_c,-temp_b(:,:,end,:,:,:));
else
    tTV_update = 0;
end

end