function PSF_visualizer(narm_frames, kspace_info, fov_factor, ifsave)
    % Visualize PSF of the spiral trajectory for RTHawk data for a set of
    % user defined number of arms per frame. Uses NUFFT.
    % 
    % Note!
    % This function is for dual TE data, and it performs splitting on the
    % kspace trajectory. If single TE data exists, splitting is not
    % necessary. A future update to this function will add dual TE as a
    % 'case' rather than a built-in functionality.
    arguments
        narm_frames = [10,5,2]
        kspace_info = []
        fov_factor = 1
        ifsave = 1
    end
    
    N = length(narm_frames);
    
    for i = 1:N
        narm_frame = narm_frames(i);
        [PSF, PSF_cropped, mainlobe_line, para] = RTHawk_spiral_PSF(narm_frame, kspace_info, fov_factor);
        
        PSF = PSF(:,:,end);
        PSF_cropped = PSF_cropped(:,:,end);
        mainlobe_line = mainlobe_line(:,end);
        
        if i == 1
            PSF_all = zeros([size(PSF), N]);
            PSF_cropped_all = zeros([size(PSF_cropped), N]);
            mainlobe_line_all = zeros([size(mainlobe_line), N]);

            PSF_all(:,:,1) = PSF;
            PSF_cropped_all(:,:,1) = PSF_cropped;
            mainlobe_line_all(:,1) = mainlobe_line;
        else
            PSF_all(:,:,i) = PSF;
            PSF_cropped_all(:,:,i) = PSF_cropped;
            mainlobe_line_all(:,i) = mainlobe_line;
        end
    end
    
    %% plot PSF image and save it
    
    % declare axes
    axis_data_high = (para.Recon.FOV * 10 / 2) * fov_factor; % upper bound for new XTick vector
    axis_data_low = - axis_data_high; % lower bound for new XTick vector
    axis_data = round(linspace(axis_data_low, axis_data_high, size(PSF, 1))); % new XTickLabel vector
    
    for i = 1:N
        % plot
        subplot(1,N,i)
        imagesc('XData', axis_data, 'YData', axis_data, 'CData', log10(abs(PSF_all(:,:,i))))
        caxis([-5.5 0])
        xlabel('X (cm)'); ylabel('Y (cm)'); axis image; colorbar; colormap gray
        title(['Log-scale PSF for ', num2str(narm_frames(i)), ' arms per frame'])
    end
    
    input('Resize the figure, then hit enter in the Command Window.')
    
    if ifsave
        saveas(gcf, './figures/logscale_PSF.png')
    end
    
    %% zoom into mainlobe
    % declare axes
    axis_data_high = (para.Recon.FOV * 100 / 2) * fov_factor * size(PSF_cropped, 1) / size(PSF, 1); % upper bound for new XTick vector
    axis_data_low = - axis_data_high; % lower bound for new XTick vector
    axis_data = round(linspace(axis_data_low, axis_data_high, size(PSF_cropped, 1))); % new XTickLabel vector
    
    for i = 1:N
        % plot
        subplot(1,N,i)
        imagesc(axis_data, axis_data, log10(abs(PSF_cropped_all(:,:,i))));
        caxis([-2 0])
        xlabel('X (mm)'); ylabel('Y (mm)'); axis image; colorbar; colormap gray
        title(['Mainlobe of PSF for ', num2str(narm_frames(i)), ' arms per frame (log-scale)'])
    end
    
    if ifsave
        saveas(gcf, './figures/mainlobe_of_PSF.png')
    end
    
    %% plot mainlobe line profile
    close all
    for i = 1:N
        % plot
        plot(axis_data, log10(abs(mainlobe_line_all(:,i))));
        hold on
    end
    
    str = '';
    for i = 1:N
        if i == 1
            str = [str, num2str(narm_frames(i))]; %#ok<*AGROW>
        elseif i < N
            str = [str, ', ', num2str(narm_frames(i))];
        else
            str = [str, ', and ', num2str(narm_frames(i))];
        end
    end
    title(['Mainlobe of PSF for ', str, ' arms per frame (log-scale)'])
    xlabel('X (mm)'); ylabel('PSF intensity');
    
    legend_labels = cell([1 N]);
    for i = 1:N
        legend_labels{i} = num2str(narm_frames(i));
    end
    legend(legend_labels)
    
    input('Resize the figure, then hit enter in the Command Window.')
    
    if ifsave
        saveas(gcf, './figures/line_mainlobe_of_PSF.png')
    end
end