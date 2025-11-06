% Photon shot noise
R = 50; % Ohm
ee = 1.60217662e-19; % C
kb = 1.380649e-23; % J/K
T = 300;% K
resp = 0.53; % A/W
gain = 1e4; % V/W
conv_gain = resp*gain/2; %V/W /2 for 50 Ohm load
P = 0.03e-3 ;% W optical power
I = resp*P; %A
BW = 45e6; %Hz
h = 6.62607015e-34;                 %  [J / Hz]
c = 2.99792458e8;                   %  [m / s]
probe_wavelength = 795e-9;          %   [m]
hf = h * c / probe_wavelength;      %   [J]
quantum_efficiency = 0.9;


% noise V RMS
V_SN_RMS = conv_gain*(2*ee*I*BW)^0.5/resp;
V_TH_RMS = (4*kb*T*R*BW)^0.5;

% Noise PSD V/srHz
PSD_SN = V_SN_RMS/(BW^0.5);
PSD_TH = V_TH_RMS/(BW^0.5);

% Noise PSN V/srHz (roy)
PSD_SN_roy = conv_gain * sqrt(2 * hf * P) / quantum_efficiency; % needs to add quantum efficiency


% Signal
V_sig = conv_gain*P;

% Minimum NEP PDB450A
min_NEP_103 = 123e-12*0.53e3/2; % V/srHz
min_NEP_104 = 28.9e-12*0.53e4/2; % V/srHz


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% LIA quantization noise ENBW

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bits = 16;
time_constants = [100e-6, 300e-6, 1e-3, 3e-3, 10e-3, 30e-3, 100e-3...
    300e-3, 1];                         % s
enbw = 3 ./ (32 * time_constants);      % Hz
sensitivity = 1;                        % Vrms

% noise V RMS
V_QUANT_RMS = sensitivity / (2 ^ bits);

% Noise PSD V/srHz
PSD_QUAN = 2 * sqrt(2) * V_QUANT_RMS ./ sqrt(enbw);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% white noise simulation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

noise_amplitude = 2;                            % V
sinamp = 200;                                   % V
fs = 100000;                                    % samples / s
meas_time = 1;                                  % s
points = fs / meas_time;
tt = linspace(0, meas_time, points);            % s
fpass = 312;                                    % Hz
noise_power = noise_amplitude^2 * fs / 2;       % V^2 * Hz
white_noise = normrnd(0, sqrt(noise_power), [1, points]);
sig = sinamp * sin(2* pi * 129 * tt);
noisy_sig = sig + white_noise;
sig_filt = lowpass(noisy_sig, fpass, fs);
[Pxx, f] = pwelch(sig_filt, [], [], [], fs);

figure
plot(tt, sig)
hold on
plot(tt, noisy_sig)
plot(tt, white_noise)
legend('original', 'noisy', 'noise')

figure
loglog(f, sqrt(Pxx))
hold on
loglog(f, noise_amplitude * ones(size(f)))
legend('sig', 'amp')
grid
