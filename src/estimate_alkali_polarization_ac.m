% ========================================================================
%   This script estimates the absolute alkali polarization (in percents)
%   by applying an AC Bx magnetic field, measuring the amplitude of the
%   resultant probe signal for a range of pump powers, and fitting a
%   function of the form y=x/(x+b);
%   The frequency of the Bx magnetic is low enough for the polarization
%   to effectively reach a steady-state, but high enough to filter DC
%   noises.

%% To do : 1. restore system to state before measuremnt
%          2. No need to measure PD when system is locked
% ========================================================================

%% Closing open figures and clear command window

% close all;
% clc;

%% Waveform generator configuration

% logger.log('AC alkali polarization measurement initialized.')
% if ~flg.is_part_of_loop
%     disp('[+] Starting alkali polarization measurement');
% end

% Disable unnecessary outputs
AG2.OutputOFF(2); % By
AG3.OutputOFF(1); % ESR
AG1.OutputOFF(1); % NMR 129
AG2.OutputOFF(1); % NMR 131
% AG6.OutputOFF(1); % Trigger
% AG6.OutputOFF(2); % Something else

% Turn on AC Bx magnetic field
T_Bx = 8e-3; % sec
f_Bx = 1 / T_Bx; % Hz
A_Bx = 2; % Volts
AG1.OutputON(2);
AG1.Sin(2, f_Bx, A_Bx, 0, 0);

%% Scope configuration

% Define scopes
pump_scope = scope4;
probe_scope = scope2;

% Configure scopes to take data
% FIXME: make sure all channels are on (current throws exception if not)
chnl = 1;
Tscale = 5e-3; % sec/div
Vscale_pump = 0.3; % volt/div
Vscale_probe = 0.4; % volt/div
offset_probe = 0.0; % volt
offset_pump = 3*Vscale_pump; % volt
timeout = 1; % sec
tDelay = 10 * T_Bx; % sec

for scope=[pump_scope, probe_scope]
    % Coupling
    scope.setChCoupling(chnl, 'DC');
    scope.setChImpedance(chnl, 'M'); % 50 Ohm (F), or 1 MOhm (M)

    % Horizontal
    scope.setTscale(Tscale);
    scope.setTdelay(0.00 * Tscale);
    
    % Trigger
    scope.TrigMode('NORM');
end

% Vertical scale
pump_scope.setVscale(chnl, Vscale_pump);
probe_scope.setVscale(chnl, Vscale_probe);

% Vertical offfset
pump_scope.setVoffset(chnl, offset_pump);
probe_scope.setVoffset(chnl, offset_probe);

%% Measuring probe voltage for different AOM DC voltages

% Photo idoe voltage range
V_PD_MIN =0.1;% -0.25; % Volts
V_PD_MAX = 1.7;% 0.25; % Volts
V_PD = linspace(V_PD_MIN, V_PD_MAX, prm.num_voltage_samples);

% Shuffle vector to mitigate systematic error
perm_idx = randperm(length(V_PD));
V_PD = V_PD(perm_idx);
iperm_idx = 1:length(V_PD);
iperm_idx(perm_idx) = 1:length(V_PD);

Vpump_mean = zeros(1, numel(V_PD));
Vprobe_mean = zeros(1, numel(V_PD));
Vpump_std = zeros(1, numel(V_PD));
Vprobe_std = zeros(1, numel(V_PD));

% Iterate over PD voltages
t0 = tic;
for i=1:numel(V_PD)
    % FIXME: Replace with a normal timebar. Suggestion: use genericTimebar from cold\Utils\Misc
    if ~flg.is_part_of_loop
        genericTimebar(i, numel(V_PD), [], t0, [], 0) 
    end
    
    % Feed DC voltage to AOM
    AG3.DC(2, V_PD(i));
    
    % Wait a fixed time delay
    pause(tDelay);
    
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
    
    % Cosine wave model for the noisy data
    mdl = fittype(@(A,f,p,D,t) A*cos(2*pi*f*t+p)+D, 'independent', 't', 'coefficients', {'A','f','p','D'});
    
    % Initial parameter estimates
    A_init = 0.5 * (max(v) - min(v));
    f_init = f_Bx;
    p_init = 0;
    D_init = mean(v);
    inits = [A_init, f_init, p_init, D_init];
    
    % Lower bounds
    A_low = 0.2 * A_init;
    f_low = 0.9 * f_init;
    p_low = -pi;
    D_low = D_init - A_init;
    lows = [A_low, f_low, p_low, D_low];
    
    % Upper bounds
    A_up = 1.8 * A_init;
    f_up = 1.1 * f_init;
    p_up = pi;
    D_up = D_init + A_init;
    ups = [A_up, f_up, p_up, D_up];
    
    % Fit model
    cosine_fitresult = fit(t.', v.', mdl, 'StartPoint', inits, ...
                           'Lower', lows, 'Upper', ups);
    A = cosine_fitresult.A;
    f = cosine_fitresult.f;
    p = cosine_fitresult.p;
    D = cosine_fitresult.D;
    
    % Compute 1*sigma confidence intervals
    fitresult_confint = confint(cosine_fitresult, erf(1 / sqrt(2)));
    fitresult_confint(isnan(fitresult_confint)) = 0; % Handle NaNs
    fitresult_confint = diff(fitresult_confint, 1) / 2;
    
    Vprobe_mean(i) = A;
    Vprobe_std(i) = fitresult_confint(1);
    if Vprobe_std(i) == 0 % Handle failed interval estimations
        Vprobe_std(i) = std(v) / sqrt(numel(t));
    end

end

% logger.log('alkali polarization measurement done.')

%% Plotting data

% Sort data after scrambling
% Vpump_mean = Vpump_mean(iperm_idx);
% Vpump_std = Vpump_std(iperm_idx);
% Vprobe_mean = Vprobe_mean(iperm_idx);
% Vprobe_std = Vprobe_std(iperm_idx);

% Voltage-to-power calibration constants (see calibrate_AOM.m)
v2p_p1 = 85/1.7;%3.54; % mW/V
v2p_p2 = 0;%2.31; % mW

% Compute pump power from scope voltage using calibration constants
Ppump = v2p_p1 * Vpump_mean + v2p_p2; % mW

meas=create_fig(meas,3); hold on; box on; grid on;
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
    if gof.rsquare > 0.98
        fit_successful = 1;
    end
    fitresult_int = confint(fitresult, erf(1 / sqrt(2)));
    fitresult_int = diff(fitresult_int, 1) / 2;
    delta_b = fitresult_int(2);
    gof_rsquare = gof.rsquare;
    fit_cnt = fit_cnt + 1;
end

if fit_successful
    a = fitresult.a;
    b = fitresult.b;
    c = fitresult.c;

    % Compute polarization at typical pump power
    P0 = 85; % mW
    PRb = P0 ./ (P0 + b); % Percent
    delta_PRb = P0 ./ (P0 + b).^2 * delta_b;

    Ppump_fit = linspace(min(Ppump), max(Ppump), 1e3);
    Vprobe_fit = a*Ppump_fit./(Ppump_fit+b) + c;

    plot(Ppump_fit, Vprobe_fit, 'k--', 'LineWidth', 1);
    legend('Data', 'Fit', 'Location', 'NorthWest');
    lgd = legend;
    lgd.FontSize = 12;
    title([datestr(datetime), newline, ...
        sprintf('Polarization @ %d mW: %s%%', P0, numerr2str(PRb*100, delta_PRb*100))], ...
        ...sprintf('Polarization @ %d mW: %.2f \x00B1 %.2f%%', P0, PRb*100, delta_PRb*100)], ...
        'FontSize', 12);

%     gof
else
    logger.log('Unsuccessful polarization esimation. Data is too noisy.');
end

%% Saving useful data

meas.Vpump_mean    = Vpump_mean;   % Mean pump voltage
meas.Vpump_std     = Vpump_std;    % Std. dev. of pump voltage
meas.Vprobe_mean   = Vprobe_mean;  % Mean probe voltage
meas.Vprobe_std    = Vprobe_std;   % Std. dev. of probe voltage
meas.P0            = P0;           % Fitted P0 coefficient in P0/(P0+b)
meas.b             = b;            % Fitted b coefficient in P0/(P0+b)
meas.b_err         = delta_b;      % Error in b coefficient
meas.PRb           = PRb;          % Estimated alkali polarization @ 20mW
meas.PRb_err       = delta_PRb;    % Error in alkali polarization estimate

%% Resetting outputs
V_PD_RST = 1.7; % V
AG3.DC(2, V_PD_RST);
AG1.OutputOFF(1);