daq_data = meas;
V = prm;

% data
t = daq_data.time_s;                        % [s]
ga = (-7.44e7  )/(2*pi)/(1e4);              % [Hz/G] Xe129 gyromagnetic ratio
gb = (+2.2056e7)/(2*pi)/(1e4);              % [Hz/G] Xe131 gyromagnetic ratio
x = ga / gb;


% World_rotation_Hz = Xe131_modulation_amplitude / PID_scale / (1 - ratio_of_gyromagnetics) * PID_error_out
W_Hz = V.df131 / 5 / (1 - 1 / x) * daq_data.Det131_PIDb_V; % [Hz]
W_deg_hr = W_Hz * 360 * 3600;               % [deg/hr]
dBz = 588e-6 * daq_data.Bz_V;               % [G] using 0.588uG/V calibrated by Roy on 11.04.21     
Omega_deg_hr = W_deg_hr;                    % [deg/hr] 
Omega = Omega_deg_hr/360/3600;              % [Hz] 

% Plotting the measured PID outputs (world rotation and dBz)
figure(111); 
subplot(2, 2, 1);
hold on
plot(t,(dBz-mean(dBz))*1e6*1e3,'color',[0 0 0.8]);
plot(t,Omega_deg_hr-mean(Omega_deg_hr),'color',[0.6 0.2 0]);
legend('\DeltaB_z[nG]','\Omega [^\circdeg/hr]');
set(gca,'fontsize',12); grid on;
xlabel('time (sec)','fontsize',12);
ylabel('\Omega [^\circdeg/hr]','fontsize',12);


% noise analyzation theta
S(1, :) = dBz; S(2,:) = Omega;               % substitute signal
fs = 1 / mean(diff(t));                      % sampling frequency
[dBz_PSD, f] = pwelch(dBz(2:end) - mean(dBz(2:end)), [], [], [], fs); % calculate Bz correction noise-PSD % [n]=(rad/s)^2/Hz
dOmega_PSD = pwelch(Omega(2:end)- mean(Omega(2:end)), [], [], [], fs); % calculate inertial rotation noise-PSD
N(1, :) = sqrt(dBz_PSD); 
N(2, :) = sqrt(dOmega_PSD);                  % substitute noise % [N]=Hz/rtHz
N_mdeg_rthr = sqrt(dOmega_PSD) * 360 * sqrt(3600) * 1e3;  % [mdeg / rthr]



bw = 10; %   [Hz]
f_bw = f(f < bw);
N_mdeg_rthr_bw = N_mdeg_rthr(f < bw);
f_bw = f_bw(2:end);
N_mdeg_rthr_bw = N_mdeg_rthr_bw(2:end);

% fit PSD
noise_fit = fit_psd(f_bw, N_mdeg_rthr_bw);
noise_fit_conf = confint(noise_fit, erf(1 / sqrt(2)));
noise_fit_conf = diff(noise_fit_conf, 1) / 2;

% generate a piecewise model
f_fit = f_bw(2:end)';
noise_parts = ones(5, length(f_fit));
noise_parts(1, :) = noise_parts(1, :) .* noise_fit.a2 .* (f_fit.^2);
noise_parts(2, :) = noise_parts(2, :) .* noise_fit.a1 .* (f_fit);
noise_parts(3, :) = noise_parts(3, :) .* noise_fit.a0;
noise_parts(4, :) = noise_parts(4, :) .* noise_fit.am1_2 .* (1 ./ sqrt(f_fit));
noise_parts(5, :) = noise_parts(5, :) .* noise_fit.am1 .* (1 ./ f_fit);

% get piecewise of noise model fit
noise_piecewise = max(noise_parts, [], 1);
noise_fit_eval = feval(noise_fit, f_fit)';
noise_fit_piecewise = [noise_fit_eval, nan, noise_piecewise];
noise_f_fit = [f_fit, nan, f_fit];

% fit name for legend
fit_name = [num_and_err(noise_fit.a2, noise_fit_conf(1)),    '\cdotf^2 + ',...
            num_and_err(noise_fit.a1, noise_fit_conf(2)),    '\cdotf^1 + ',...
            num_and_err(noise_fit.a0, noise_fit_conf(3)),    '\cdotf^0 + ', ...
            num_and_err(noise_fit.am1_2, noise_fit_conf(4)), '\cdotf^{-0.5} + ', ...
            num_and_err(noise_fit.am1, noise_fit_conf(5)),   '\cdotf^{-1}'];



try
    [sigy, tau] = (ADEV(Omega_deg_hr(2:end), fs)); % calculate inertial rotation Allan Var % [sigy]=(deg/hr), [tau]=s
catch; tau = nan; sigy = nan;
end

% plot noise-PSD
figure(111); 
subplot(2, 2, [3, 4]);
% loglog(f, N(1, :) * 1e6, f, N_mdeg_rthr);
loglog(f, N_mdeg_rthr, 'color', [0 0.4470 0.7410]);
hold on

% plot fit and spline
loglog(f_fit, noise_fit_eval, ':', 'linew', 2, 'color', [0 0.2470 0.7410]);
loglog(f_fit, noise_piecewise, 'linew', 2, 'color', [0 0.2470 0.7410]);


set(gca, 'fontsize', 12);
xlabel('frequency (Hz)', 'fontsize', 12); 
ylabel('\delta\omega (mdeg/\surdhr)', 'fontsize', 12); 
grid on;
ARW = 0.005 * 116e-6; %  sqrt(2 * pi) / (360* 60) = 116e-6
AWN = 0.012 / 360; 
BiasIns = 0.02 / 360 / 3600; % [deg/hr]/360/3600 = [Hz^2 / Hz]
fH = 0.5; % could be up to 5...
sqrth2 = 2 * pi * AWN / sqrt(3 * fH); 
sqrth0 = sqrt(2) * ARW; 
sqrthM1 = BiasIns / sqrt(2 * log(2));
nPSD_NGC = sqrt((sqrth2 * f).^2 + (sqrth0 * ones(size(f))).^2 + (sqrthM1./sqrt(f)).^2 );
hold all; 
plot(f, nPSD_NGC * (1e3 / 116e-6), ':k', 'linew', 3);
legend('\Omega', fit_name, 'piecewise', 'NGC');
axis tight

% adev
subplot(2, 2, 2);
loglog(tau, sigy, 'linew', 2);
grid on;
xlabel('\tau (s)');
ylabel('\sigma_{\ity} (deg/hr)'); 
title('Allan deviation');
set(gca, 'XMinorGrid', 'off')
set(gca, 'yMinorGrid', 'off')

[ind, ff] = min(find(tau > 1));
ARW = sigy(ind) ./ sqrt(tau);
AWN = sigy(ind) ./ tau;
hold on;
plot(tau, ARW, ':r', 'linewidth', 2)
plot(tau, AWN, ':g', 'linewidth', 2)
ADAV_NGC = sqrt((0.012 ./ (tau / 3600)).^2 + (0.005 ./ sqrt(tau / 3600)).^2 + 0.02^2);
plot(tau, ADAV_NGC, '--k', 'linew', 2);
legend('Measured','ARW','AWN','NGC');
tau = 1;







