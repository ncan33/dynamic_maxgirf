function writeSpinBenchWaveforms(wfFile, waveforms, includeGradient, includeRF)
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

numAxes = 0;
if(includeGradient)
    numAxes = numAxes + 3;
end
if(includeRF)
    numAxes = numAxes + 2;
end
% if(numAxes == 0) 
%     disp("Data file must have either gradient or RF!");
%     return;
% end


% if(fid == -1) 
%     disp(msg);
%     return;
% end

% change this to 'ieee-le' if using little-endian
% waveforms = reshape(waveforms, [length(waveforms)/numAxes numAxes]);

fid = fopen(wfFile, 'w');
fwrite(fid, waveforms, 'int16', 'ieee-be');
fclose(fid);
% waveforms = reshape(alldata, [length(alldata)/numAxes numAxes]);
