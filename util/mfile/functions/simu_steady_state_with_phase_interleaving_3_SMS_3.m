function SI_one_slice = simu_steady_state_with_phase_interleaving_3_SMS_3(flip_angle_slice,phase_slice,Mz,M0,TR,Slab_T1,Nalpha)
% Nalpha = 1:3000; % number of alpha pulses

set = repmat([1,2,3],[1,Nalpha/3]);
phase = vec(repmat([0,1,2],[3,Nalpha/9]))';
SI_set = zeros(3,Nalpha/3);


for i=1:Nalpha
    flip_angle_all = zeros(2600,1);
    
    SMS_slices = [1:1000,601:1600,1201:2200];
    if set(i) == 2
        SMS_slices = SMS_slices + 200;
    elseif set(i) == 3
        SMS_slices = SMS_slices + 400;
    end
    
    flip_angle_all(SMS_slices(1:1000)) = flip_angle_all(SMS_slices(1:1000)) + flip_angle_slice;
    flip_angle_all(SMS_slices(1001:2000)) = flip_angle_all(SMS_slices(1001:2000)) + flip_angle_slice;
    flip_angle_all(SMS_slices(2001:3000)) = flip_angle_all(SMS_slices(2001:3000)) + flip_angle_slice;

    Mz(:,end+1) = Mz(:,end).*cos(flip_angle_all/180*pi);
    
    switch phase(i)
        case 0
            n_phase = [0,0,0];
        case 1
            n_phase = [0,1,2];
        case 2
            n_phase = [0,2,1];
    end
    
    for iii = [1,4,7]-1+set(i)
        n_SMS_slice = ceil(iii/3);
        N = ceil(i/3);
        slice_location = (1:1000) + 200*(iii-1);

        SI_set(set(i),N) = SI_set(set(i),N) + sum(Mz(slice_location,end-1).*sin(flip_angle_slice/180*pi).*phase_slice) .*exp(1i*n_phase(n_SMS_slice)*2*pi/3);
        
    end

%     for iii=1:9
%         slice_location = (1:1000) + 200*(iii-1);
%         SI_one_slice(iii,i) = sum(Mz(slice_location,end-1).*sin(flip_angle_slice/180*pi).*phase_slice);
%     end
    Mz(:,end+1) = M0 - (M0 - Mz(:,end)).*exp(-TR./Slab_T1);

end

SI_one_slice(1,:) = abs(SI_set(1,1:3:end) + SI_set(1,2:3:end) + SI_set(1,3:3:end))/3;
SI_one_slice(2,:) = abs(SI_set(2,1:3:end) + SI_set(2,2:3:end) + SI_set(2,3:3:end))/3;
SI_one_slice(3,:) = abs(SI_set(3,1:3:end) + SI_set(3,2:3:end) + SI_set(3,3:3:end))/3;

SI_one_slice(4,:) = abs(SI_set(1,1:3:end)+SI_set(1,2:3:end)*exp(-1i*2*pi/3)+SI_set(1,3:3:end)*exp(-1i*4*pi/3))/3;
SI_one_slice(5,:) = abs(SI_set(2,1:3:end)+SI_set(2,2:3:end)*exp(-1i*2*pi/3)+SI_set(2,3:3:end)*exp(-1i*4*pi/3))/3;
SI_one_slice(6,:) = abs(SI_set(3,1:3:end)+SI_set(3,2:3:end)*exp(-1i*2*pi/3)+SI_set(3,3:3:end)*exp(-1i*4*pi/3))/3;

SI_one_slice(7,:) = abs(SI_set(1,1:3:end)+SI_set(1,2:3:end)*exp(-1i*4*pi/3)+SI_set(1,3:3:end)*exp(-1i*2*pi/3))/3;
SI_one_slice(8,:) = abs(SI_set(2,1:3:end)+SI_set(2,2:3:end)*exp(-1i*4*pi/3)+SI_set(2,3:3:end)*exp(-1i*2*pi/3))/3;
SI_one_slice(9,:) = abs(SI_set(3,1:3:end)+SI_set(3,2:3:end)*exp(-1i*4*pi/3)+SI_set(3,3:3:end)*exp(-1i*2*pi/3))/3;

% 
% SI_one_slice_temp([1,4,7],:) = SI_one_slice([1,4,7],1:3:end);
% SI_one_slice_temp([2,5,8],:) = SI_one_slice([2,5,8],2:3:end);
% SI_one_slice_temp([3,6,9],:) = SI_one_slice([3,6,9],3:3:end);

% SI_one_slice = SI_one_slice_temp;
end