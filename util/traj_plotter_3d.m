%% Plot trajectory in 3D
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

%% View order
view_order = kspace_info.viewOrder;
new_view_order = view_order;
counter = 0;
for i = 1:size(view_order, 2)
    if view_order(i) > 10
        new_view_order(i) = view_order(i) - 10;
        counter = counter + 1;
    end
end

plot(new_view_order(1:50))
title('First 50 elements of the view order')
