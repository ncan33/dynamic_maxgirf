function [kx_echo_1, ky_echo_1, kx_echo_2, ky_echo_2] = traj_plotter(kspace_info, narm_frame)
% note: if you scroll, you will see the dual TE processed kspace trajectory

    arguments
        kspace_info
        narm_frame = 10
    end

    %% process traj based on view_order
    nview = kspace_info.user_interleaves; % total number of interleaves in original kspace trajectory
    view_order = kspace_info.viewOrder;

    echo_idx = view_order > nview;
    echo_idx = echo_idx + 1; % if the element in echo_idx == 1, then echo_1,
                             % and when echo_idx == 2, then echo_2

    nsample = kspace_info.extent(1); % number of samples
    ncoil = kspace_info.extent(2); % number of coils

    narm_total = min(sum(echo_idx == 1), sum(echo_idx == 2));
    nframes = floor(narm_total / narm_frame);
    narm_total = nframes * narm_frame;

    view_order_echo_1 = view_order(echo_idx == 1); % view_order for echo_1
    view_order_echo_2 = view_order(echo_idx == 2); % view_order for echo_2
    view_order_echo_1(narm_total + 1 : end) = []; % discard excess views
    view_order_echo_2(narm_total + 1 : end) = []; % discard excess views
    view_order_echo_1  = reshape(view_order_echo_1, [narm_frame, nframes]);
    view_order_echo_2  = reshape(view_order_echo_2, [narm_frame, nframes]);

    kx = kspace_info.kx_GIRF; % kx
    ky = kspace_info.ky_GIRF; % ky

    kx_echo_1 = zeros(nsample, narm_frame, nframes); %kx_echo_1 = zeros(size(kspace_echo_1,1), size(kspace_echo_1,2), size(kspace_echo_1,3));
    kx_echo_2 = kx_echo_1;
    ky_echo_1 = kx_echo_1;
    ky_echo_2 = kx_echo_1;
    
    for i = 1:narm_frame
        for j = 1:nframes
            kx_echo_1(:,i,j) = kx(:, view_order_echo_1(i,j));
            kx_echo_2(:,i,j) = kx(:, view_order_echo_2(i,j) - 10);

            ky_echo_1(:,i,j) = ky(:, view_order_echo_1(i,j));
            ky_echo_2(:,i,j) = ky(:, view_order_echo_2(i,j) - 10);
        end
    end
    
    %% plot
    clf
    hold on
    axis square
    xlim([-0.5 0.5])
    ylim([-0.5 0.5])
    
    %{
    for frame = 1:2
        for arm = 1:10
            plot(kx_echo_2(:,arm,frame), ky_echo_2(:,arm,frame));
            %xlim([-0.30, -0.10])
            %ylim([0.10 0.30])
            title(['k-space trajectory, arm ', num2str(arm),' and frame ', num2str(frame)]);
            pause(1);
        end
        cla
        hold on
    end
    %}
    
    counter = 0;
    for frame = 1:15
        for arm = 1:10
            counter = counter + 1;
            plot(kx_echo_2(:,arm,frame), ky_echo_2(:,arm,frame));
            %xlim([-0.30, -0.10])
            %ylim([0.10 0.30])
            title(['k-space trajectory for frame ', num2str(counter)]);
            pause(0.1);
            cla
            hold on
        end
    end
end