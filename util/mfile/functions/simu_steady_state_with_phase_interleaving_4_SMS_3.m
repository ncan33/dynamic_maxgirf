function SI_one_slice = simu_steady_state_with_phase_interleaving_4_SMS_3(flip_angle_slice,phase_slice,Mz,M0,TR,Slab_T1,Nalpha)
% Nalpha = 1:3000; % number of alpha pulses

set = repmat([1,2,3,4],[1,Nalpha/4]);
phase = vec(repmat([0,1,2],[4,Nalpha/12]))';
SI_set = zeros(4,Nalpha/4);


for i=1:Nalpha
    flip_angle_all = zeros(3200,1);
    
    SMS_slices = [1:1000,801:1800,1601:2600];
    if set(i) == 2
        SMS_slices = SMS_slices + 200;
    elseif set(i) == 3
        SMS_slices = SMS_slices + 400;
    elseif set(i) == 4
        SMS_slices = SMS_slices + 600;
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
    
    for iii = [1,5,9]-1+set(i)
        n_SMS_slice = ceil(iii/4);
        N = ceil(i/4);
        slice_location = (1:1000) + 200*(iii-1);

        SI_set(set(i),N) = SI_set(set(i),N) + sum(Mz(slice_location,end-1).*sin(flip_angle_slice/180*pi).*phase_slice) .*exp(1i*n_phase(n_SMS_slice)*2*pi/3);
        
    end

%     for iii=1:9
%         slice_location = (1:1000) + 200*(iii-1);
%         SI_one_slice(iii,i) = sum(Mz(slice_location,end-1).*sin(flip_angle_slice/180*pi).*phase_slice);
%     end
    Mz(:,end+1) = M0 - (M0 - Mz(:,end)).*exp(-TR./Slab_T1);

end

SI_one_slice(1,:) = abs(SI_set(1,1:4:end) + SI_set(1,2:4:end) + SI_set(1,3:4:end))/3;
SI_one_slice(2,:) = abs(SI_set(2,1:4:end) + SI_set(2,2:4:end) + SI_set(2,3:4:end))/3;
SI_one_slice(3,:) = abs(SI_set(3,1:4:end) + SI_set(3,2:4:end) + SI_set(3,3:4:end))/3;
SI_one_slice(4,:) = abs(SI_set(4,1:4:end) + SI_set(4,2:4:end) + SI_set(4,3:4:end))/3;

SI_one_slice(5,:) = abs(SI_set(1,1:4:end)+SI_set(1,2:4:end)*exp(-1i*2*pi/3)+SI_set(1,3:4:end)*exp(-1i*4*pi/3))/3;
SI_one_slice(6,:) = abs(SI_set(2,1:4:end)+SI_set(2,2:4:end)*exp(-1i*2*pi/3)+SI_set(2,3:4:end)*exp(-1i*4*pi/3))/3;
SI_one_slice(7,:) = abs(SI_set(3,1:4:end)+SI_set(3,2:4:end)*exp(-1i*2*pi/3)+SI_set(3,3:4:end)*exp(-1i*4*pi/3))/3;
SI_one_slice(8,:) = abs(SI_set(4,1:4:end)+SI_set(4,2:4:end)*exp(-1i*2*pi/3)+SI_set(4,3:4:end)*exp(-1i*4*pi/3))/3;

SI_one_slice(9,:) = abs(SI_set(1,1:4:end)+SI_set(1,2:4:end)*exp(-1i*4*pi/3)+SI_set(1,3:4:end)*exp(-1i*2*pi/3))/3;
SI_one_slice(10,:) = abs(SI_set(2,1:4:end)+SI_set(2,2:4:end)*exp(-1i*4*pi/3)+SI_set(2,3:4:end)*exp(-1i*2*pi/3))/3;
SI_one_slice(11,:) = abs(SI_set(3,1:4:end)+SI_set(3,2:4:end)*exp(-1i*4*pi/3)+SI_set(3,3:4:end)*exp(-1i*2*pi/3))/3;
SI_one_slice(12,:) = abs(SI_set(4,1:4:end)+SI_set(4,2:4:end)*exp(-1i*4*pi/3)+SI_set(4,3:4:end)*exp(-1i*2*pi/3))/3;


end