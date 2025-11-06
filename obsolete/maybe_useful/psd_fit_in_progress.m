% [https://en.wikipedia.org/wiki/Allan_variance#Measurement_issues]

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   We are linearly fitting the loglog srPSD data with y = a*x + b for  
%   different a's in order to extract the next physical quantities: 
%   AWN (angle white noise):        a = 1
%   ARW (angle random walk):        a = 0
%   RBI (rate bias instability):    a = -0.5
%
%
%   Why is this the right think to do ???
%   Following the Allan Variance wikipedia page, the PSD fit in units of    
%   deg^2/hr^2/Hz in the power coefficients notation for AWN, ARW and RBI 
%   is given by
%
%                       S = h2*f^2 + h0 + h-1*f^-1
%
%   The above single fit can be broken down into three linear fits in
%   loglog presentation of the sqrt PSD in different parts of the spectrum 
%   because different parts of the spectrum are dominated by different
%   types of noises, e.g. the awn fit would be
%
%                       S = h2*f^2 + h0 + h-1*f^-1
%
%   ==>         log10(sqrt(S_awn)) = log10(sqrt(h2*f^2))
%               log10(sqrt(S_awn)) = 0.5*log10(h2*f^2)
%               log10(sqrt(S_awn)) = 0.5*log10(h2) + 0.5*log10(f^2)
%               log10(sqrt(S_awn)) = 0.5*log10(h2) + log10(f)
%
%   ==>         y = log10(sqrt(S_awn))
%               x = log10(f)
%               a = 1
%               b = 0.5*log10(h2)
%
%   when       h2 = 2*(2*pi*PHI)^2
%
%   where the factor 2 comes from the fact that we are fitting the single
%   sided PSD instead of the PSD itself.
%
%   ==>         PHI(IEEE) = 10^(b) / (2*pi) / sqrt(2)
%
%
%   
%   Next, we compute these quantities from the offset of the fit (b) as 
%   follows:
%               PHI = 10^psd_awn_fit.b / (2*pi) / sqrt(2)
%               N   = 10^psd_awn_fit.b / sqrt(2)
%               B   = 10^psd_rbi_fit.b * sqrt(2*pi) / sqrt(2)
%
%   In the case where the sqrt(PSD) has the units of mdeg/hr/srHz we get
%   the following units for the noise quantities above
%
%               [PHI] = mdeg/hr/Hz^(3/2)
%               [N]   = mdeg/hr/srHz
%               [B]   = mdeg/hr
%
%   Finally, we would like to transform the units to Walker's notations 
%   in order to compare our results with the NGC's results
%
%        [PHI_walker] = deg
%        [N_walker]   = deg/srhr
%        [B_walker]   = deg/hr
%
%   Let us put aside for a moment the AWN (PHI) transformation adn focus on
%   the other two quantaties transformation that would be as follows
%
%        N_walker = N / 60  <== 1/srHz = sqrt(s) = sqrt(h/3600) = srhr/60
%        B_walker = B
%
%   The problem with PHI is that it cannot be transformed into deg units.
%   We deal with that by talking about AWN (PHI) but normalizing the units 
%   as if we are fitting AQN (angle quantization noise) which is denoted
%   in the IEEE notation as Q. The AQN term in the PSD fit looks like this
%
%                  AQN_term = 4Q^2 / tau * sin(pi*f*tau)^2
% 
%   In the general case tau would be my instruments bandwidth (Nyquist)
%   while f is equal to where I am measuring in the spectrum, e.g. in our 
%   system (after smoothing the data) tau = 1/fs = 0.01s, where f~1Hz so 
%   we are working in a regime where tau * f << 1 therefore
%
%                  AQN_term ~ 2*tau*(2*pi*f*Q)^2    :(single sided PSD)
%
%   the units of the AQN are
%
%                               [Q] = deg/hr/Hz
%
%   Notice that the units of Q can be transformed into deg as follow
%
%                          Q_walker = Q / 3600
%
%   It is important to note that in the regime of our system both AQN and 
%   AWN have the same slope in the PSD which is f^2. Therefore, we cannot
%   distiguish between the two quantatites in the data. Nevertheless, we
%   work with the AWN quantity and the units shuold be carfully exameaned.
%   
%          
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parameters
bw = 1; %   [Hz]
fH = 0.5; % could be in [0.5, 5] 

% Constants
ga = (-7.44e7  )/(2*pi)/(1e4);              % [Hz/G] Xe129 gyromagnetic ratio
gb = (+2.2056e7)/(2*pi)/(1e4);              % [Hz/G] Xe131 gyromagnetic ratio
x = ga / gb;

% data
daq_data = meas;
V = prm;
t = daq_data.time_s;                        % [s]


% World_rotation_Hz = Xe131_modulation_amplitude / PID_scale / (1 - ratio_of_gyromagnetics) * PID_error_out
W_Hz = V.df131 / 5 / (1 - 1 / x) * daq_data.Det131_PIDb_V;      % [Hz]
Omega_deg_hr = W_Hz * 360 * 3600;                               % [deg/hr]
dBz = 588e-6 * daq_data.Bz_V;                                   % [G] using 0.588uG/V calibrated by Roy on 11.04.21     


% noise analyzation theta
S = zeros(2, length(dBz));
S(1, :) = dBz; S(2,:) = Omega_deg_hr;               % substitute signal
fs = 1 / mean(diff(t));                             % sampling frequency
[dBz_PSD, f] = pwelch(dBz(2:end) - mean(dBz(2:end)), [], [], [], fs);                               % calculate Bz correction noise-PSD % [G^2/Hz]
dOmega_PSD_deg2_hr2_hz = pwelch(Omega_deg_hr(2:end)- mean(Omega_deg_hr(2:end)), [], [], [], fs);    % calculate inertial rotation noise-PSD [deg^2/hr^2/Hz]
dOmega_PSD_deg2_hr = dOmega_PSD_deg2_hr2_hz / 3600; % [deg^2/hr] --> (hr = 3600s = 3600/Hz)
N = zeros(2, length(dBz_PSD));
N(1, :) = sqrt(dBz_PSD);                            % [G/rtHz]
N(2, :) = sqrt(dOmega_PSD_deg2_hr);                 % [deg/rthr]
N_mdeg_rthr = sqrt(dOmega_PSD_deg2_hr) * 1e3;       % [mdeg/rthr]



%% Analizing only data in Bandwidth

f_bw = f(f < bw);
f_bw = f_bw(2:end);
dOmega_PSD_deg2_hr_bw = dOmega_PSD_deg2_hr(f < bw);     % [deg^2/hr]
dOmega_PSD_deg2_hr_bw = dOmega_PSD_deg2_hr_bw(2:end);   % [deg^2/hr]


%% fit PSD

% setting up the frequency limits
psd_rbi_lim_s = 1e-2;                    % [Hz]
psd_rbi_lim_e = 10 * psd_rbi_lim_s;      % [Hz]
psd_arw_lim_s = 1e-1;                    % [Hz]
psd_arw_lim_e = 3 * psd_arw_lim_s;       % [Hz]
psd_awn_lim_s = bw / 4;                  % [Hz]
psd_awn_lim_e = bw;                      % [Hz]

dOmega_PSD_mdeg_hr_srHz_log = log10(sqrt(dOmega_PSD_deg2_hr2_hz) * 1e3);
f_log = log10(f);

% uplim = max(dOmega_PSD_mdeg_hr_srHz_log);
% lowlim = min(dOmega_PSD_mdeg_hr_srHz_log);
% figure(123)
% plot(f_log, dOmega_PSD_mdeg_hr_srHz_log)
% hold on
% plot(log10([psd_rbi_lim_s psd_rbi_lim_s]), [lowlim uplim], '--k', 'LineWidth', 1)
% plot(log10([psd_rbi_lim_e psd_rbi_lim_e]), [lowlim uplim], '--k', 'LineWidth', 1)
% plot(log10([psd_arw_lim_s psd_arw_lim_s]), [lowlim uplim], '--k', 'LineWidth', 1)
% plot(log10([psd_arw_lim_e psd_arw_lim_e]), [lowlim uplim], '--k', 'LineWidth', 1)
% plot(log10([psd_awn_lim_s psd_awn_lim_s]), [lowlim uplim], '--k', 'LineWidth', 1)
% plot(log10([psd_awn_lim_e psd_awn_lim_e]), [lowlim uplim], '--k', 'LineWidth', 1)

% verify that f, and PSD are both row vectors 
f_size = size(f);
psd_size = size(dOmega_PSD_mdeg_hr_srHz_log);
if f_size(1) > 1
    f = f';
end
if psd_size(1) > 1
    dOmega_PSD_mdeg_hr_srHz_log = dOmega_PSD_mdeg_hr_srHz_log';
end



% Data division and Fit in loglog
% AWN fit
psd_awn_bool = all([f > psd_awn_lim_s; f < psd_awn_lim_e]);
f_awn_linear = f(psd_awn_bool);
f_awn_log = log10(f_awn_linear);
dOmega_PSD_mdeg_hr_srHz_log_awn = dOmega_PSD_mdeg_hr_srHz_log(psd_awn_bool);
psd_awn_fit_op = fitoptions('Method','NonlinearLeastSquares',...
                        'StartPoint', (dOmega_PSD_mdeg_hr_srHz_log_awn(1) - 1*f_awn_log(1)));
psd_awn_noise_model = fittype('1 * x + b', ...
                          'coefficients', {'b'},...
                          'options', psd_awn_fit_op);
psd_awn_fit = fit(f_awn_log', dOmega_PSD_mdeg_hr_srHz_log_awn', psd_awn_noise_model);
psd_awn_fit_conf = confint(psd_awn_fit, erf(1 / sqrt(2)));
psd_awn_fit_conf = diff(psd_awn_fit_conf, 1) / 2;

% ARW fit
psd_arw_bool = all([f > psd_arw_lim_s; f < psd_arw_lim_e]);
f_arw_linear = f(psd_arw_bool);
f_arw_log = log10(f_arw_linear);
dOmega_PSD_mdeg_hr_srHz_log_arw = dOmega_PSD_mdeg_hr_srHz_log(psd_arw_bool);
psd_arw_fit_op = fitoptions('Method','NonlinearLeastSquares',...
                        'StartPoint', (dOmega_PSD_mdeg_hr_srHz_log_arw(1)));
psd_arw_noise_model = fittype('0 * x + b', ...
                          'coefficients', {'b'},...
                          'options', psd_arw_fit_op);
psd_arw_fit = fit(f_arw_log', dOmega_PSD_mdeg_hr_srHz_log_arw', psd_arw_noise_model);
psd_arw_fit_conf = confint(psd_arw_fit, erf(1 / sqrt(2)));
psd_arw_fit_conf = diff(psd_arw_fit_conf, 1) / 2;


% RBI fit
psd_rbi_bool = all([f > psd_rbi_lim_s; f < psd_rbi_lim_e]);
f_rbi_linear = f(psd_rbi_bool);
f_rbi_log = log10(f_rbi_linear);
dOmega_PSD_mdeg_hr_srHz_log_rbi = dOmega_PSD_mdeg_hr_srHz_log(psd_rbi_bool);
psd_rbi_fit_op = fitoptions('Method','NonlinearLeastSquares',...
                        'StartPoint', (dOmega_PSD_mdeg_hr_srHz_log_rbi(1) + 0.5 * f_rbi_log(1)));
psd_rbi_noise_model = fittype('-0.5 * x + b', ...
                          'coefficients', {'b'},...
                          'options', psd_rbi_fit_op);
psd_rbi_fit = fit(f_rbi_log', dOmega_PSD_mdeg_hr_srHz_log_rbi', psd_rbi_noise_model);
psd_rbi_fit_conf = confint(psd_rbi_fit, erf(1 / sqrt(2)));
psd_rbi_fit_conf = diff(psd_rbi_fit_conf, 1) / 2;


% Generating linear fit data for plot
psd_awn_fit_plot = 10.^feval(psd_awn_fit, f_awn_log) / 60;      % [mdeg/srhr]
psd_arw_fit_plot = 10.^feval(psd_arw_fit, f_arw_log) / 60;      % [mdeg/srhr]
psd_rbi_fit_plot = 10.^feval(psd_rbi_fit, f_rbi_log) / 60;      % [mdeg/srhr]


% Extract physical quantities from fit (IEEE notation)
psd_PHI_mdeg_hr_Hz3_2 = 10^psd_awn_fit.b / 2 / pi / sqrt(2);    % [mdeg/hr/Hz^(3/2)]
psd_N_mdeg_hr_srHz = 10^psd_arw_fit.b / sqrt(2);                % [mdeg/hr/srHz]
psd_B_mdeg_hr = 10^psd_rbi_fit.b * sqrt(2 * pi) / sqrt(2);      % [mdeg/hr]

% Extract physical quantities from fit (Walker's notation)
psd_PHI_deg_srHz = psd_PHI_mdeg_hr_Hz3_2 / 3600 * 1e-3;         % [deg/srHz] (in Walker, the units of this term is deg)
psd_N_deg_srhr = psd_N_mdeg_hr_srHz / 60 * 1e-3;                % [deg/srhr]
psd_B_deg_hr = psd_B_mdeg_hr * 1e-3;                            % [deg/hr]

% Extract errors from fit (Walker notation) (needs to be double checked, I did this in a rush)
psd_PHI_err = (10^(psd_awn_fit.b + psd_awn_fit_conf) - 10^(psd_awn_fit.b)) / 2 / pi / sqrt(2) / 3600 * 1e-3;
psd_N_err = (10^(psd_arw_fit.b + psd_arw_fit_conf) - 10^(arw_allan_fit.b)) / sqrt(2) / 60 * 1e-3;
psd_B_err = (10^(psd_rbi_fit.b + psd_rbi_fit_conf) - 10^(rbi_allan_fit.b)) * sqrt(2 * pi) / sqrt(2) * 1e-3;


% % Plot fit
% plot(f_awn_log, feval(psd_awn_fit, f_awn_log)', '--m')
% plot(f_arw_log, feval(psd_arw_fit, f_arw_log)', '--b')
% plot(f_rbi_log, feval(psd_rbi_fit, f_rbi_log)', '--r')


% w = 2 * pi * f;
% gen_psd_mdeg_hr_srHz = psd_PHI_mdeg_hr_Hz .* w + psd_N_mdeg_hr_srHz .* ones(size(w)) + psd_B_mdeg_hr ./ (w.^0.5);  % [mdeg/hr/srHz]
% gen_psd_mdeg_srhr = gen_psd_mdeg_hr_srHz / 60;  % [mdeg/srhr]



 
 

%% Fit Allen

% compute the Allan Variance of signal
[avar, tau1] = (AVAR(Omega_deg_hr(2:end), fs));

% cleaning the datat from NaNs
avar_bool = ~isnan(avar);
avar_clean = avar(avar_bool);   % [deg^2/h^2]
tau1 = tau1(avar_bool);         % [s]

% Set limits for fit
awn_lim = 0.1;  % [s]
arw_lim_s = 0.5;  % [s]
arw_lim_e = 1.2;  % [s]
rbi_lim = 1.2;  % [s]

% % set nnn = 1000 to generate more samples of avar
% nnn = 1;
% a = 300;
% b = 1e3;
% avar_clean = [avar_clean a.*randn(1, nnn)+b];      % [deg^2/hr^2]
% tau1 = [tau1 logspace(1, 1, nnn)];

% work with log-log adev
avar_deg2_sec2 = avar_clean / (3600^2);             % [deg^2/s^2]
adev_log = log10(sqrt(avar_deg2_sec2));             % [deg/s]
tau2_log = log10(tau1);

% uplim = max(adev_log);
% lowlim = min(adev_log);


% Data division and Fit
% AWN fit
awn_bool = tau2_log < awn_lim;
tau_awn = tau2_log(awn_bool);
adev_awn = adev_log(awn_bool);
awn_fit_op = fitoptions('Method','NonlinearLeastSquares',...
                        'StartPoint', (adev_awn(1) + 2*tau_awn(1)));
awn_noise_model = fittype('-1 * x + b', ...
                          'coefficients', {'b'},...
                          'options', awn_fit_op);
awn_allan_fit = fit(tau_awn', adev_awn', awn_noise_model);
awn_allan_fit_conf = confint(awn_allan_fit, erf(1 / sqrt(2)));
awn_allan_fit_conf = diff(awn_allan_fit_conf, 1) / 2;

% ARW fit
arw_bool = all([tau2_log > arw_lim_s ; tau2_log < arw_lim_e]);
tau_arw = tau2_log(arw_bool);
adev_arw = adev_log(arw_bool);
arw_fit_op = fitoptions('Method','NonlinearLeastSquares',...
                        'StartPoint', (adev_arw(1) + 1*tau_arw(1)));
arw_noise_model = fittype('-0.5 * x + b', ...
                          'coefficients', {'b'},...
                          'options', arw_fit_op);
arw_allan_fit = fit(tau_arw', adev_arw', arw_noise_model);
arw_allan_fit_conf = confint(arw_allan_fit, erf(1 / sqrt(2)));
arw_allan_fit_conf = diff(arw_allan_fit_conf, 1) / 2;

% RBI fit
rbi_bool = tau2_log > rbi_lim;
tau_rbi = tau2_log(rbi_bool);
adev_rbi = adev_log(rbi_bool);
rbi_fit_op = fitoptions('Method','NonlinearLeastSquares',...
                        'StartPoint', (adev_rbi(1)));
rbi_noise_model = fittype('-0 * x + b', ...
                          'coefficients', {'b'},...
                          'options', rbi_fit_op);
rbi_allan_fit = fit(tau_rbi', adev_rbi', rbi_noise_model);
rbi_allan_fit_conf = confint(rbi_allan_fit, erf(1 / sqrt(2)));
rbi_allan_fit_conf = diff(rbi_allan_fit_conf, 1) / 2;

% Extract physical quantities from fit (IEEE notation)
allan_PHI_deg_srHz = 10^awn_allan_fit.b / sqrt(3) / sqrt(fH);       % [deg/srHz]
allan_N_deg_s_srHz = 10^arw_allan_fit.b;                            % [deg/s/srHz] = [deg/srs]
allan_B_deg_s = 10^rbi_allan_fit.b * sqrt(pi / (2 * log(2)));       % [deg/s]

% Extract physical quantities for piecewise
allan_awn_deg = 10^awn_allan_fit.b;                                 % [deg]         PHI*sr(fH)
allan_arw_deg_s_srHz = 10^arw_allan_fit.b;                          % [deg/s/srHz]
allan_rbi_deg_s = 10^rbi_allan_fit.b;                               % [deg/s]

% Extract physical quantities from fit (Walker notation)
allan_PHI_deg_srHz = 10^awn_allan_fit.b / sqrt(3) / sqrt(fH);           % [deg/srHz] (in Walker, the units of this term is deg)
allan_N_deg_srhr = 10^arw_allan_fit.b * 60;                             % [deg/srhr]
allan_B_deg_hr = 10^rbi_allan_fit.b * sqrt(pi / (2 * log(2))) * 3600;   % [deg/hr]

% Extract errors from fit (Walker notation) (needs to be double checked, I did this in a rush)
allan_PHI_err = (10^(awn_allan_fit.b + awn_allan_fit_conf) - 10^(awn_allan_fit.b)) / sqrt(3) / sqrt(fH);
allan_N_err = (10^(arw_allan_fit.b + arw_allan_fit_conf) - 10^(arw_allan_fit.b)) * 60;
allan_B_err = (10^(rbi_allan_fit.b + rbi_allan_fit_conf) - 10^(rbi_allan_fit.b)) * sqrt(pi / (2 * log(2))) * 3600;

% Generating a Piecewise model
allen_dev_model_piecewise_parts = ones(3, length(tau1));
allen_dev_model_piecewise_parts(1, :) = allen_dev_model_piecewise_parts(1, :) .* allan_awn_deg ./ tau1;                     % [deg/s]
allen_dev_model_piecewise_parts(2, :) = allen_dev_model_piecewise_parts(2, :) .* allan_arw_deg_s_srHz ./ (tau1.^(0.5));     % [deg/s]
allen_dev_model_piecewise_parts(3, :) = allen_dev_model_piecewise_parts(3, :) .* allan_rbi_deg_s;                           % [deg/s]


% Getting the piecewise of allan model fit
adev_deg_h_piecewise = max(allen_dev_model_piecewise_parts, [], 1) * 3600;      % [deg/h]



%% 
% calculate inertial rotation Allan Var % [sigy]=(deg/hr), [tau]=s
try 
    [sigy, tau] = (ADEV(Omega_deg_hr(2:end), fs)); 
catch; tau = nan; sigy = nan;
end

% NGC specs (Spin-Exchange Pumped NMR Gyros - arXiv:1604.03982)
ngc_awn = 0.012;    % [deg]
ngc_arw = 0.005;    % [deg/srhr]
ngc_bias = 0.02;    % [deg/hr]

%% True up until here

% Normalizing to world rotation pure units [Hz^2/Hz]
deg_srhr2rad_hz = sqrt(2 * pi) / 360 / (sqrt(3600));
ngc_h2 = (2 * pi * ngc_awn / 360 / sqrt(3 * fH))^2;             % [1/Hz] 
ngc_h0 = (sqrt(2) * ngc_arw * sqrt(2 * pi) / 360 / sqrt(3600))^2;               % [Hz^2/Hz]
ngc_hm1 = (ngc_bias / sqrt(2 * log(2)) / 360 / 3600)^2;         % [Hz^2/Hz]

% Generating NGC spec PSD line in sqrt([Hz^2/Hz]) units
ngc_psd = sqrt(ngc_h2 * f.^2 + ngc_h0 * f.^0  + ngc_hm1 * f.^(-1));
ngc_psd_mdeg_srhr = ngc_psd * 1e3 / sqrt(2 * pi) * 360 * sqrt(3600);

% Allan Deviation of NGC spec
ADAV_NGC = sqrt((ngc_awn ./ (tau / 3600)).^2 + (ngc_arw ./ sqrt(tau / 3600)).^2 + ngc_bias^2);
% ADAV_gen = sqrt((PHI_deg_srhr_Hz ./ (tau / 3600)).^2 + (N_deg_srhr ./ sqrt(tau / 3600)).^2 + B_deg_rthr_rtHz^2); % / sqrt(360)

%% Setting up the legend text

awn_text = ['AWN (\Phi)  [deg/srHz]', ' | ' 'Allan: ', num2str(allan_PHI_deg_srHz), ' | ', 'Psd: ', num2str(psd_PHI_deg_srHz), ' | '];
arw_text = ['ARW (N)  [deg/srhr]', '   | ' 'Allan: ', num2str(allan_N_deg_srhr), ' | ', 'Psd: ', num2str(psd_N_deg_srhr), ' | '];
rbi_text = ['RBI   (B)   [deg/hr]  ', '   | ' 'Allan: ', num2str(allan_B_deg_hr), ' | ', 'Psd: ', num2str(psd_B_deg_hr), ' | '];


%% Plotting 
figure1 = figure(113);

% Plot measured data
axes1 = axes('Parent', figure1, 'Position', [0.13 0.65 0.33 0.27]);
hold(axes1,'on');
plot(t,(dBz-mean(dBz))*1e6*1e3,'color',[0 0 0.8], 'Parent', axes1);
plot(t,Omega_deg_hr-mean(Omega_deg_hr),'color',[0.6 0.2 0], 'Parent', axes1);
legend('\DeltaB_z[nG]','\Omega [^\circdeg/hr]');
grid(axes1, 'on');
set(axes1, 'FontSize',10);
xlabel('time (sec)','fontsize',10);
ylabel('\Omega [^\circdeg/hr]','fontsize',10);

% Plot PSD
axes2 = axes('Parent', figure1, 'Position', [0.13 0.11 0.60 0.37]);
hold(axes2,'on');
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
xlabel('frequency (Hz)', 'fontsize', 12); 
ylabel('\delta\omega (mdeg/\surdhr)', 'fontsize', 12); 
loglog(f, N_mdeg_rthr, 'color', [0 0.4470 0.7410], 'Parent', axes2);
plot(f, ngc_psd_mdeg_srhr, ':k', 'linew', 3, 'Parent', axes2);
%loglog(f, gen_psd_mdeg_srhr, 'Parent', axes2);

% Plot PSD fit
plot(f_awn_linear, psd_awn_fit_plot, ':g', 'linew', 2, 'Parent', axes2)
plot(f_arw_linear, psd_arw_fit_plot, ':r', 'linew', 2, 'Parent', axes2)
plot(f_rbi_linear, psd_rbi_fit_plot, ':m', 'linew', 2, 'Parent', axes2)

box(axes2,'on');
grid(axes2,'on');
axis(axes2,'tight');
legend2 = legend('\Omega', 'NGC', awn_text, arw_text, rbi_text);
set(legend2, 'Position', [0.8 0.11 0.08 0.37]);


% Plot Allan Deviation
axes3 = axes('Parent', figure1, 'Position', [0.57 0.65 0.33 0.27]);
loglog(tau1, sqrt(avar_clean), 'linew', 2, 'Parent', axes3);
hold(axes3,'on');
plot(tau, ADAV_NGC, '--k', 'linew', 2, 'Parent', axes3)
plot(tau1(awn_bool), allen_dev_model_piecewise_parts(1, awn_bool) * 3600, ':g', 'linew', 2, 'Parent', axes3)
plot(tau1(arw_bool), allen_dev_model_piecewise_parts(2, arw_bool) * 3600, ':r', 'linew', 2, 'Parent', axes3)
plot(tau1(rbi_bool), allen_dev_model_piecewise_parts(3, rbi_bool) * 3600, ':m', 'linew', 2, 'Parent', axes3)

% aaaaa = ngc_awn * ones(size(tau1)) * 3600;
% plot(tau1(awn_bool), aaaaa(awn_bool)./tau1(awn_bool), ':r', 'linew', 2, 'Parent', axes3)

% plot(10.^[awn_lim awn_lim], 10.^[lowlim uplim] * 3600, '--k', 'LineWidth', 1)
% plot(10.^[arw_lim_s arw_lim_s], 10.^[lowlim uplim] * 3600, '--k', 'LineWidth', 1)
% plot(10.^[arw_lim_e arw_lim_e], 10.^[lowlim uplim] * 3600, '--k', 'LineWidth', 1)
% plot(10.^[rbi_lim rbi_lim], 10.^[lowlim uplim] * 3600, '--k', 'LineWidth', 1)
%legend('Measured', 'NGC', 'AWN', 'ARW', 'RBI');
grid on;
xlabel('\tau (s)', 'fontsize',10);
ylabel('\sigma_{\ity} (deg/hr)', 'fontsize',10); 
title('Allan deviation');
set(gca, 'XMinorGrid', 'off')
set(gca, 'yMinorGrid', 'off')


