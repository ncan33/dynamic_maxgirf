kx = kspace_info.kx;
ky = kspace_info.ky;
readout_time = kspace_info.user_readoutTime;
time_vector = linspace(0, readout_time, size(kx, 1));


plot3(kx(:,1), time_vector, ky(:,1))
xlabel('kx')
ylabel('Time (ms)')
zlabel('ky')
title('k-space trajectory (1 spiral arm only)')

figure
plot3(kx(:,1:5), repmat(time_vector.', 1, 5), ky(:,1:5))
xlabel('kx')
ylabel('Time (ms)')
zlabel('ky')
title('k-space trajectory (5 spiral arms)')

figure
plot3(kx, repmat(time_vector.', 1, 10), ky)
xlabel('kx')
ylabel('Time (ms)')
zlabel('ky')
title('k-space trajectory (all 10 spiral arms)')