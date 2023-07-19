function [time, Gx, Gy, Gz, Rf, RfPhase] = loadRthWaveForm(file_name, ifplot)

fid = fopen(file_name);
i = 1;
while ~feof(fid)
    temp = fgetl(fid);
    if temp(1) ~= '%'
        WaveForm(:, i) = sscanf(temp, '%f %f %f %f %f %f');
        i = i + 1;
    end
end
fclose(fid);

time = WaveForm(1, :);
Gx = WaveForm(2, :);
Gy = WaveForm(3, :);
Gz = WaveForm(4, :);
Rf = WaveForm(5, :);
RfPhase = WaveForm(6, :);

if ifplot
    figure
    subplot(4, 1, 1)
    plot(time, Rf)
    subplot(4, 1, 2)
    plot(time, Gx)
    subplot(4, 1, 3)
    plot(time, Gy)
    subplot(4, 1, 4)
    plot(time, Gz)
    
    figure
    subplot(4, 1, 1)
    plot(Rf)
    subplot(4, 1, 2)
    plot(Gx)
    subplot(4, 1, 3)
    plot(Gy)
    subplot(4, 1, 4)
    plot(Gz)
end
    
    