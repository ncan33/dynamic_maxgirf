function SI_one_slice = simu_steady_state_interleaving_3_SMS_3(flip_angle_slice,Mz,M0,TR,Slab_T1,Nalpha)
% Nalpha = 1:3000; % number of alpha pulses

interleaving_slices = repmat([1,2,3],[1,Nalpha/3]);


for i=1:Nalpha
    flip_angle_all = zeros(2600,1);
    
    SMS_slices = [1:1000,601:1600,1201:2200];
    if interleaving_slices(i) == 2
        SMS_slices = SMS_slices + 200;
    elseif interleaving_slices(i) == 3
        SMS_slices = SMS_slices + 400;
    end
    
    flip_angle_all(SMS_slices(1:1000)) = flip_angle_all(SMS_slices(1:1000)) + flip_angle_slice;
    flip_angle_all(SMS_slices(1001:2000)) = flip_angle_all(SMS_slices(1001:2000)) + flip_angle_slice;
    flip_angle_all(SMS_slices(2001:3000)) = flip_angle_all(SMS_slices(2001:3000)) + flip_angle_slice;

    Mz(:,end+1) = Mz(:,end).*cos(flip_angle_all/180*pi);

    for iii=1:9
        slice_location = (1:1000) + 200*(iii-1);
        SI_one_slice(iii,i) = sum(Mz(slice_location,end-1).*sin(flip_angle_slice/180*pi));
    end
    Mz(:,end+1) = M0 - (M0 - Mz(:,end)).*exp(-TR./Slab_T1);

end

SI_one_slice_temp([1,4,7],:) = SI_one_slice([1,4,7],1:3:end);
SI_one_slice_temp([2,5,8],:) = SI_one_slice([2,5,8],2:3:end);
SI_one_slice_temp([3,6,9],:) = SI_one_slice([3,6,9],3:3:end);

SI_one_slice = SI_one_slice_temp;
end