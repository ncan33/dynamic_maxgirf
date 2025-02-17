function [kPred, GPred, PhOut] = apply_GIRF_yt(gradients_nominal, dt, R, acqDelay, gradDelay)
% function [kPred, GPred] = apply_GIRF(gradients_nominal, dt, R, tRR)
%
% % UNITS
% gradients_nominal [ G/cm ]
% dt                [ s ]
%
% Hack to R to handle field strength (backwards compatibility)
% R.R = rotation matrix;
% R.T = field strength {}
%
% acqDelay in us: this should be matched to the acqusition delay applied during the acquisition
% negative value means ADC preceeds by the absolute amount of acq_delay
% postive value means ADC starts the absolute amount of acq_delay after gradient applied

field_T = 0.55;
R = R.R;
%% LOAD GIRF (Field/scanner dependent)
% Field selection needs to be handled upstream, or include the headers
% in to this file

% girf_file = '/Volumes/GoogleDrive/My Drive/3.USC/Research/Codes/rthawk/GIRF/GIRF_20210827.mat';
% girf_file = '/Volumes/GoogleDrive/My Drive/3.USC/Research/Codes/rthawk/GIRF/GIRF_20210829.mat';
girf_file = 'GIRF_20210829.mat';

% === Load file ===
try
    girf_data = load(girf_file);
%     disp(['Using ' girf_file]);
    girf0 = girf_data.girf.girf0;
    girf1 = girf_data.girf.girf1;
  
    girf.freq = girf_data.girf.freq;
    dtGIRF = 1/(girf_data.girf.freq(end)-girf_data.girf.freq(1));
    dtGIRF = 10e-6;
catch
    error ('Couldn''t find the GIRF file');
end

dtSim = dt;
l_GIRF = length(girf1);
[samples, interleaves, gs] = size(gradients_nominal);

%%   GIRF prediction   %
if nargin < 4
    acqDelay = 0;
end

%%
clear G0 GNom GPred kNom kPred

% allocation
Nominal     = zeros([samples, interleaves, 3], 'double');
Predicted   = zeros([samples, interleaves, 3], 'double');
GNom        = zeros(samples, 3, interleaves, 'double');
GPred       = zeros(samples, 3, interleaves, 'double');
kNom        = zeros(samples, 3, interleaves, 'double');
kPred       = zeros(samples, 3, interleaves, 'double');

B0Pred      = zeros([samples, interleaves, 3], 'double');
PhPred      = zeros([samples, interleaves, 3], 'double');
PhOut       = zeros(samples, interleaves, 'double');

% GIRF process
G0 = gradients_nominal;

%Rotate into physical coordinates
G0 = reshape(G0, [samples * interleaves, gs]);
G0 = (R * G0.').';
G0 = reshape(G0, [samples, interleaves, gs]);

gradDelayInPhysics = gradDelay.';
%--Loop through x,y,z gradient trajectories--%
for ax = 1:gs
    ADCshift = acqDelay*1e-6 + (-4e-6);
    % if we want to use the estimated gradDelay from RTHawk
    % ADCshift = -6*1e-6 - gradDelayInPhysics(ax)*1e-6 ;
    
    % Zeropad in time domain to match frequency resolution of GIRF (match readout length)
    L = round(dtGIRF * l_GIRF / dtSim); % when waveform not at GRT
    G = zeros([L, interleaves], 'double');
    
    G(1:samples, :) = G0(:, :, ax);
    
    % Make a waveform periodic by returning to zero
    H = G(samples, :) .* hanning(400);
    G(samples + (1:200), :) = H(201:end, :);
    
    %FFT nominal gradient
    dw = 1 / (L * dtSim); % frequency resolution [Hz]
    w = (-floor(L/2):ceil(L/2)-1).' * dw; % [Hz]
    
    GIRF1 = interp1(girf.freq,girf1(:,ax),w);
    GIRF1(isnan(GIRF1))=0;
    
    GIRF0 = interp1(girf.freq,girf0(:,ax),w);
    GIRF0(isnan(GIRF0))=0;
    
    I = fftshift(fft(ifftshift(G, 1)), 1);
    P = I .* GIRF1 .* exp(1j * ADCshift * 2 * pi * w);
    
    P_b0 = I .* GIRF0 .* exp(1j * ADCshift * 2 * pi * w);
    b0_ec = real(fftshift(ifft(ifftshift(P_b0, 1)), 1) * L / L);
    
    B0Pred(:, :, ax) = b0_ec(1:samples, :);
    PhPred(:, :, ax) = 2 * pi * dt * cumsum(B0Pred(:, :, ax),1);
    
    PredGrad = fftshift(ifft(ifftshift(P, 1)), 1) * L / L;
    NomGrad  = fftshift(ifft(ifftshift(I, 1)), 1)  * L / L;
    
    Nominal(:, :, ax)   = NomGrad(1:samples, :);
    Predicted(:, :, ax) = PredGrad(1:samples, :);
end


PhOut = sum(PhPred, 3);
PhOut = PhOut * 100; % girf0 is in rad/(mT/m) whereas gradient input is in 1/100 mT/m

%rotate back to logical coordinates
Nominal     = reshape(Nominal, [samples * interleaves, gs]);
GNom        = (R.'*Nominal.').';
GNom        = reshape(GNom, [samples, interleaves, gs]);

Predicted   = reshape(Predicted, [samples * interleaves, gs]);
GPred       = (R.'*Predicted.').';
GPred       = reshape(GPred, [samples, interleaves, gs]);

GNom    = real(GNom);
GPred   = real(GPred);

%Integrate to get k-space trajectory from gradient
kNom  = cumsum(GNom);
kPred = cumsum(GPred);


% Scale k-space in units of 1/cm
gamma = 2.67522212e+8; % gyromagnetic ratio for 1H [rad/sec/T]
% Modified by NGL: 2.675e8 => gamma
kPred = 0.01*(gamma/(2*pi))*(kPred*0.01)*dt; % (kPred*0.01):: assuming gradients in are in G/cm!!!
kNom  = 0.01*(gamma/(2*pi))*(kNom*0.01)*dt; % (kPred*0.01):: assuming gradients in are in G/cm!!!

%Permute
% kPred = permute(kPred,[1 3 2]);
% kNom  = permute(kNom,[1 3 2]);
% GPred = permute(GPred,[1 3 2]);
% GNom  = permute(GNom,[1 3 2]);

% figure,
% subplot(2,2,1); plot(real(kNom(:,1,1)),'or');
% hold on; plot(real(kPred(:,1,1)),'ob');
% subplot(2,2,2); plot(real(kNom(:,1,2)),'or');
% hold on; plot(real(kPred(:,1,2)),'ob');
% 
% subplot(2,2,3); plot(real(GNom(:,1,1)),'or');
% hold on; plot(real(GPred(:,1,1)),'ob');
% subplot(2,2,4); plot(real(GNom(:,1,2)),'or');
% hold on; plot(real(GPred(:,1,2)),'ob');

end
