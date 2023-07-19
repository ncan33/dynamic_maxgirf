function [g_slice_selective] = dzgrad_single_slice(g_platue, n_points_platue, g_dwell, rf_dwell, g_max, slew_max, ifprewinder, ifplot)
%--------------------------------------------------------------------------
%   dzgrad_blipped_caipi(g_platue, n_points_platue, g_dwell, rf_dwell, g_max, slew_max)
%--------------------------------------------------------------------------
%   Design slice selective gradient waveform for bilpped-caipi
%--------------------------------------------------------------------------
%   Inputs:
%       - g_platue          [1x1] [G/cm]
%       - n_points_platue
%       - g_dwell           [1x1] [ms]
%       - rf_dwell          [1x1] [ms]
%       - g_max             [1x1] [G/cm]
%       - slew_max          [1x1] [G/cm/ms]

%--------------------------------------------------------------------------
%               |                                      |
%               |       slice selective gradient       |
%               |                                      |
%             A |        ____________                  |
%               |       /            \                 |
%             B |  ..../..............\..............  |
%               |                      \        /      |
%             C |                       \______/       |
%               |      a b           c d e     f g     |
%--------------------------------------------------------------------------
%   a-b:    time reverse of 3-4
%   b-c:    g_platue, n_points_platue
%   c-d:    g_rewinder_platue_to_0, 
%   d-e:    g_rewinder_0_to_bottom, n_step_0_to_rewinder_bottom
%   e-f:    , n_steps_rewinder_bottom
%   f-g:    time reverse of 4-5
%   A:      g_platue
%   B:      zero
%   C:      g_rewinder_bottom
%% design gradient waveform
dwell_scale = g_dwell / rf_dwell;

g_platue = ones(n_points_platue / dwell_scale, 1) * g_platue(1);
area_under_g_platue = sum(g_platue * g_dwell); % [G/cm] * [ms]

% n_steps_to_platue = g_platue(1) / slew_max / g_dwell; % [G/cm] / [G/cm/ms] / [ms] -> unitless
g_rewinder_platue_to_0 = g_platue(1) : - slew_max * g_dwell : 0;
g_rewinder_platue_to_0 = g_rewinder_platue_to_0(2:end)';
area_remining = area_under_g_platue / 2 + sum(g_rewinder_platue_to_0) * g_dwell; % [G/cm] * [ms]

g_rewinder_bottom = - g_max; 
if g_rewinder_bottom^2 / slew_max > area_remining % [G/cm] * [G/cm] / [G/cm/ms] -> [G/cm] * [ms] > [G/cm] * [ms] 
    g_rewinder_bottom = sqrt(area_remining * slew_max / 1000);
    % need more work here
else
    n_step_0_to_rewinder_bottom = - (g_rewinder_bottom - g_rewinder_platue_to_0(end)) / slew_max / g_dwell; % [G/cm] / [G/cm/ms] / [ms] -> unitless
    g_rewinder_0_to_bottom = g_rewinder_platue_to_0(end) : (g_rewinder_bottom - g_rewinder_platue_to_0(end)) / n_step_0_to_rewinder_bottom : g_rewinder_bottom;
    g_rewinder_0_to_bottom = g_rewinder_0_to_bottom(2:end)';
    area_rewinder_ramps = sum(g_rewinder_0_to_bottom) * g_dwell * 2;
    area_rewinder_platue = area_remining + area_rewinder_ramps;
    n_steps_rewinder_bottom = area_rewinder_platue / g_dwell / - g_rewinder_bottom;
    n_steps_rewinder_bottom = ceil(n_steps_rewinder_bottom);
    g_rewinder_bottom = - area_rewinder_platue / g_dwell / n_steps_rewinder_bottom;
    
    flag = g_rewinder_bottom(1) > g_rewinder_0_to_bottom(end);
    while flag 
        n_steps_rewinder_bottom = n_steps_rewinder_bottom + 2;
        area_rewinder_platue = area_rewinder_platue - g_rewinder_0_to_bottom(end) * g_dwell * 2;
        g_rewinder_0_to_bottom(end) = [];
        g_rewinder_bottom = - area_rewinder_platue / g_dwell / n_steps_rewinder_bottom;
        flag = g_rewinder_bottom(1) > g_rewinder_0_to_bottom(end);
    end
    g_rewinder_platue = ones(n_steps_rewinder_bottom, 1) * g_rewinder_bottom;
end

g_rewinder = [g_rewinder_platue_to_0; g_rewinder_0_to_bottom; g_rewinder_platue; g_rewinder_0_to_bottom(end:-1:1)];

%% slice selective gradient, at dewll time of g_dwell
g_slice_selective = [g_platue; g_rewinder];
if ifprewinder
    g_prewinder = g_rewinder(end:-1:1);
    g_slice_selective = [g_prewinder; g_slice_selective];
end


% n_points_rewinder = length(g_rewinder);

% dp = linspace(dp(1), dp(2), length(g_slice_selective_rf_dwell));

n_points_g = length(g_slice_selective);
t = g_dwell * (1:n_points_g);

%% plot
if ifplot
    figure
    plot(t, g_slice_selective, 'Color', colors(1), 'LineWidth', 2)
    set(gca, 'XLim', [t(1), t(end)])
    xlabel 'ms'
    ylabel 'G/cm'
    grid on
    set(gca, 'FontSize', 14)
end


