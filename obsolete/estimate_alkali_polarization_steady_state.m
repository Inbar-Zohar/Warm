% ========================================================================
%   This script estimates the absolute Alkali polarization (in percents)
%   by measuring the steady-state probe DC voltage signal for a range of
%   pump powers, and fitting a function of the form y=x/(x+b);
%   During the experiment both Bz and Bx magnetic fields are on.
% ========================================================================

%% Closing open figures and clear command window

% close all;
% clc;

%% Waveform generator configuration

logger.log('Steady-state Alkali polarization measurement initialized.')
disp('[+] Starting Alkali polarization measurement');

% Disable unnecessary outputs
AG2.OutputOFF(2); % By
AG3.OutputOFF(1); % ESR
AG1.OutputOFF(1); % NMR 129
AG2.OutputOFF(1); % NMR 131
AG6.OutputOFF(1); % Trigger
AG6.OutputOFF(2); % Something else

% Turn on Bx magnetic field
AG1.OutputON(2);
AG1.DC(2, 1);

%% Scope configuration

% Define scopes
pump_scope = scope4;
probe_scope = scope2;

% Configure scopes to take data
chnl = 1;
Tscale = 20e-6; % sec/div
Vscale_pump = 5; % volt/div
Vscale_probe = 0.1; % volt/div
timeout = 1; % sec
tDelay = 2e-2; % sec

for scope=[pump_scope, probe_scope]
    % Coupling
    scope.setChCoupling(chnl, 'DC');
    scope.setChImpedance(chnl, 'M'); % 50 Ohm (F), or 1 MOhm (M)

    % Horizontal
    scope.setTscale(Tscale);
    scope.setTdelay(0.00 * Tscale);

    % Vertical
    scope.setVoffset(chnl, 0);
    
    % Trigger
    scope.TrigMode('NORM');
end

% Vertical scale
pump_scope.setVscale(chnl, Vscale_pump);
probe_scope.setVscale(chnl, Vscale_probe);

%% Measuring probe voltage for different AOM DC voltages

% AOM voltage range
V_AOM = -0.25:0.01:0.25;

% Shuffle vector to mitigate systematic error
perm_idx = randperm(length(V_AOM));
V_AOM = V_AOM(perm_idx);
iperm_idx = 1:length(V_AOM);
iperm_idx(perm_idx) = 1:length(V_AOM);

Vpump_mean = zeros(1, numel(V_AOM));
Vprobe_mean = zeros(1, numel(V_AOM));
Vpump_std = zeros(1, numel(V_AOM));
Vprobe_std = zeros(1, numel(V_AOM));

% Iterate over AOM voltages
for i=1:numel(V_AOM)
    disp(i);

    % Feed DC voltage to AOM
    AG3.DC(2, V_AOM(i));
    
    % Read data from pump scope
    pump_scope.Single;
    pump_scope.ForceTrig;
    pump_scope.readyToRead(timeout);
    [t, v, ~] = pump_scope.Read(chnl);
    Vpump_mean(i) = mean(v);
    Vpump_std(i) = std(v) / sqrt(numel(t));
    
    % Read data from probe scope
    probe_scope.Single;
    probe_scope.ForceTrig;
    probe_scope.readyToRead(timeout);
    [t, v, ~] = probe_scope.Read(chnl);
    Vprobe_mean(i) = mean(v);
    Vprobe_std(i) = std(v) / sqrt(numel(t));
    
    % Wait a fixed time delay
    pause(tDelay);
end

logger.log('Alkali polarization measurement done.')

%% Plotting data

% Sort data after scrambling
Vpump_mean = Vpump_mean(iperm_idx);
Vpump_std = Vpump_std(iperm_idx);
Vprobe_mean = Vprobe_mean(iperm_idx);
Vprobe_std = Vprobe_std(iperm_idx);

% Voltage-to-power calibration constants (see calibrate_AOM.m)
v2p_p1 = 3.54; % mW/V
v2p_p2 = 2.31; % mW

% Compute pump power from scope voltage using calibration constants
Ppump = v2p_p1 * Vpump_mean + v2p_p2; % mW

figure; hold on; box on;
errorbar(Ppump, Vprobe_mean, Vprobe_std, Vprobe_std, v2p_p1*Vpump_std, v2p_p1*Vpump_std, '.', 'LineWidth', 1);
xlabel('Pump power (mW)', 'FontSize', 12);
ylabel('Probe voltage (V)', 'FontSize', 12);

%% Analyzing data

gof_rsquare = 0;
fit_successful = 0;
fit_cnt = 0;
while ~fit_successful && fit_cnt < 30 % Repeat until successful
    [fitresult, gof] = fit(Ppump.', Vprobe_mean.', 'a*x/(x+b)+c', ...
        'StartPoint', [0, 0, min(Vprobe_mean)]);
    if gof.rsquare > 0.99
        fit_successful = 1;
    end
    fitresult_int = confint(fitresult, erf(1 / sqrt(2)));
    fitresult_int = diff(fitresult_int, 1) / 2;
    delta_b = fitresult_int(2);
    gof_rsquare = gof.rsquare
    fit_cnt = fit_cnt + 1;
end

if fit_successful
    a = fitresult.a;
    b = fitresult.b;
    c = fitresult.c;

    % Compute polarization at typical pump power
    P0 = 20; % mW
    PRb = P0 ./ (P0 + b); % Percent
    delta_PRb = P0 ./ (P0 + b).^2 * delta_b;

    Ppump_fit = linspace(min(Ppump), max(Ppump), 1e3);
    Vprobe_fit = a*Ppump_fit./(Ppump_fit+b) + c;

    plot(Ppump_fit, Vprobe_fit, 'k--', 'LineWidth', 1);
    ylim([-0.3,0.15]);
    legend('Data', 'Fit', 'Location', 'NorthWest');
    lgd = legend;
    lgd.FontSize = 12;
    title([datestr(datetime), newline, ...
           sprintf('Polarization @ %d mW: %.2f \x00B1 %.2f%%', P0, PRb*100, delta_PRb*100)], ...
           'FontSize', 12);

    gof
else
    disp('[+] Unsuccessful polarization esimation. Data is too noisy.');
end

%% Resetting AOM

V_AOM_RST = 0.5; % V
AG3.DC(2, V_AOM_RST);