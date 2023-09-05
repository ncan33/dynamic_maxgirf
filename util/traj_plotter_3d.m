kx = kspace_info.kx;
ky = kspace_info.ky;
readout_time = kspace_info.user_readoutTime;
time_vector = linspace(0, readout_time, size(kx, 1));

plot3(kx(:,1), time_vector, ky(:,1))

figure

plot3(kx, repmat(time_vector.', 1, 10), ky)
