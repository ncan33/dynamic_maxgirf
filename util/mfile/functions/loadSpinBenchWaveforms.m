function waveforms = loadSpinBenchWaveforms(wfFile, includeGradient, includeRF)
% waveforms = loadSpinBenchWaveforms(wfFile, includeGradient=true, includeRF=true)
%
% This is a simple MATLAB/Octave function that loads waveform files
% exported from SpinBench.
%
% Parameters are:
%  wfFile         : filename of exported file
%  includeGradient: boolean indicating whether gradient data is
%                   included in file.  Default TRUE
%  includeRF      : boolean indicating whether RF data is
%                   included in file.  Default TRUE
%

waveforms=[];

numAxes = 0;
if(includeGradient)
    numAxes = numAxes + 3;
end
if(includeRF)
    numAxes = numAxes + 2;
end
if(numAxes == 0) 
    disp("Data file must have either gradient or RF!");
    return;
end

[fid, msg] = fopen(wfFile, 'r', 'b');
if(fid == -1) 
    disp(msg);
    return;
end

% change this to 'ieee-le' if using little-endian
alldata = fread(fid, Inf, 'int16', 'ieee-be');
fclose(fid);
waveforms = reshape(alldata, [length(alldata)/numAxes numAxes]);
