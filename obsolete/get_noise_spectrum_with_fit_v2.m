function [f, N, t, S, tau, computed_quantities] = get_noise_spectrum_with_fit_v2(daq_data, V, fit_bw)
% This function plot the measured data from the PID output with its PSD and
% ADEV. In addition it fit the PSD and ADEV and extract the AWN, ARW and
% RBI from the measured data.
% [https://en.wikipedia.org/wiki/Allan_variance#Measurement_issues]

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   We are linearly fitting the loglog srPSD data with y = a*x + b for  
%   different a's in order to extract the next physical quantities: 
%   AWN (angle white noise):        a = 1
%   ARW (angle random walk):        a = 0
%   RBI (rate bias instability):    a = -0.5
%
%
%   Why is this the right thing to do ???
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
%   work with the AWN quantity (in psd), however, in ADEV we are extracting 
%   the AQN (or AWN with fH=1) so the units shuold be carfully exameaned.
%   
%          
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Analysis
switch nargin
    case 2
        fit_bw = 1;     % [Hz] - the system's bandwidth  
end

if ~isfield(V,'american_lock')
    V.american_lock = 0;
end
               
% Constants
bw = fit_bw;
dTheta_dV = 18/1;% deg/Volt
fH_ngc = 1;   % could be in [0.5, 5] - see fH in wiki Allan Variance
ga = (-7.44e7  )/(2*pi)/(1e4);              % [Hz/G] Xe129 gyromagnetic ratio
gb = (+2.2056e7)/(2*pi)/(1e4);              % [Hz/G] Xe131 gyromagnetic ratio
x = ga / gb;

% data
t = daq_data.time_s;                        % [s]

if V.american_lock
    % american lock from improvement of stability
    warning('american lock now')
    R = abs(x);
W_Hz_131 = V.df131 / 5 * daq_data.Det131_PIDb_V;      % [Hz]
W_Hz_129 = V.df131 / 5  * daq_data.Bz_V;      % [Hz]
W_Hz = W_Hz_129/(1+R)-W_Hz_131*R/(1+R);
Omega_deg_hr = W_Hz * 360 * 3600;  
dBz = (W_Hz_129+W_Hz_131)/(abs(ga)+gb);                               % [G] 100uA/V*40uA/150uA*3G/A! 

% Temp calculation for unlocked system 
% Gamma129 = 1/3.3;
% Gamma131 = 1/11;
% theta_deg_131 = dTheta_dV*daq_data.Xe131_theta_V; % Xe 131 Phase [deg]
% theta_deg_129 = dTheta_dV*daq_data.Xe129_theta_V; % Xe 131 Phase [deg]
% Omega_deg_hr = (Gamma129*theta_deg_129/(1+R)-Gamma131*theta_deg_131*R/(1+R))*3600; %  
% dBz = (Gamma129*theta_deg_129+Gamma131*theta_deg_131)*3600;%/(abs(ga)+gb);  % [G] 100uA/V*40uA/150uA*3G/A! 

else
    % World_rotation_Hz = Xe131_modulation_amplitude / PID_scale / (1 - ratio_of_gyromagnetics) * PID_error_out
W_Hz = V.df131 / 5 / (1 - 1 / x) * daq_data.Det131_PIDb_V;      % [Hz]
Omega_deg_hr = W_Hz * 360 * 3600;                               % [deg/hr]
dBz = -100e-6*40/150*3 * daq_data.Bz_V;                               % [G] 100uA/V*40uA/150uA*3G/A! 

    
end
% noise analyzation theta
S = zeros(2, length(dBz));
S(1, :) = dBz; S(2,:) = Omega_deg_hr;               % substitute signal
fs = 1 / mean(diff(t));                             % sampling frequency
%!!!! Why are we taking (2:end) and not the whole data?  !!!!!!!!!!!!!!!!!
[dBz_PSD, f] = pwelch(dBz(2:end) - mean(dBz(2:end)), [], [], [], fs);                               % calculate Bz correction noise-PSD % [G^2/Hz]
dOmega_PSD_deg2_hr2_hz = pwelch(Omega_deg_hr(2:end)- mean(Omega_deg_hr(2:end)), [], [], [], fs);    % calculate inertial rotation noise-PSD [deg^2/hr^2/Hz]
dOmega_PSD_deg2_hr = dOmega_PSD_deg2_hr2_hz / 3600; % [deg^2/hr] --> (hr = 3600s = 3600/Hz)
N = zeros(2, length(dBz_PSD));
N(1, :) = sqrt(dBz_PSD);                            % [G/rtHz]
N(2, :) = sqrt(dOmega_PSD_deg2_hr);                 % [deg/rthr]
N_mdeg_rthr = sqrt(dOmega_PSD_deg2_hr) * 1e3;       % [mdeg/rthr]

% Collect interesting calculated values
psd_PHI_deg_srHz = 0;
psd_N_deg_srhr = 0;
psd_B_deg_hr = 0;
allan_PHI_deg_srHz = 0;
allan_N_deg_srhr = 0;
allan_B_deg_hr = 0;

% from PSD 
computed_quantities.psd_PHI_deg_srHz = 0;
computed_quantities.psd_N_deg_srhr = 0;
computed_quantities.psd_B_deg_hr = 0;
computed_quantities.psd_PHI_deg_srHz_err = 0;
computed_quantities.psd_N_deg_srhr_err = 0;
computed_quantities.psd_B_deg_hr_err = 0;

% from Allan 
computed_quantities.allan_PHI_deg_srHz = 0;
computed_quantities.allan_N_deg_srhr = 0;
computed_quantities.allan_B_deg_hr = 0;
computed_quantities.allan_PHI_deg_srHz_err = 0;
computed_quantities.allan_N_deg_srhr_err = 0;
computed_quantities.allan_B_deg_hr_err = 0;

%% fit PSD

dOmega_PSD_mdeg_hr_srHz_log = log10(sqrt(dOmega_PSD_deg2_hr2_hz) * 1e3);
% taking the minimal psd value for arw window calculation
[~, psd_argmin] = min(dOmega_PSD_mdeg_hr_srHz_log);

% if the minimal psd value is at the end of spectrum, compute arw
% window according to fit_bw parameter
if psd_argmin >= length(f) / 2
   [~, psd_argmin] = min(abs(fit_bw / 10 - f));
end

% setting up the frequency limits
psd_arw_lim_e = f(psd_argmin) * 1.5;        % [Hz]
psd_arw_lim_s = f(psd_argmin) / 1.5;        % [Hz]
psd_awn_lim_s = psd_arw_lim_e;              % [Hz]
psd_awn_lim_e = fit_bw;                     % [Hz]
psd_rbi_lim_s = 1e-2;                       % [Hz]
psd_rbi_lim_e = psd_arw_lim_s;              % [Hz]



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

if sum(psd_awn_bool) > 1
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
    
    % Generating linear fit data for plot
    psd_awn_fit_plot = 10.^feval(psd_awn_fit, f_awn_log) / 60;      % [mdeg/srhr]
    
    % Extract physical quantities from fit (IEEE notation)
    psd_PHI_mdeg_hr_Hz3_2 = 10^psd_awn_fit.b / 2 / pi / sqrt(2);    % [mdeg/hr/Hz^(3/2)]
    
    % Extract physical quantities from fit (Walker's notation)
    psd_PHI_deg_srHz = psd_PHI_mdeg_hr_Hz3_2 / 3600 * 1e-3 * (2 * pi) * sqrt(2 * pi)/10;   % [deg/srHz] (in Walker, the units of this term is deg)

    % Extract errors from fit (Walker notation) (needs to be double checked, I did this in a rush)
    psd_PHI_err = (10^(psd_awn_fit.b + psd_awn_fit_conf) - 10^(psd_awn_fit.b)) / 2 / pi / sqrt(2) / 3600 * 1e-3 * (2 * pi) * sqrt(2 * pi)/10;
        warning('added /10 because I think that the extraction has error')
    % Save quantities
    computed_quantities.psd_PHI_deg_srHz = psd_PHI_deg_srHz;
    computed_quantities.psd_PHI_deg_srHz_err = psd_PHI_err;
end

% ARW fit
psd_arw_bool = all([f > psd_arw_lim_s; f < psd_arw_lim_e]);

if sum(psd_arw_bool) > 1
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
    
    % Generating linear fit data for plot
    psd_arw_fit_plot = 10.^feval(psd_arw_fit, f_arw_log) / 60;      % [mdeg/srhr]
    
    % Extract physical quantities from fit (IEEE notation)
    psd_N_mdeg_hr_srHz = 10^psd_arw_fit.b / sqrt(2);                % [mdeg/hr/srHz]

    % Extract physical quantities from fit (Walker's notation)
    psd_N_deg_srhr = psd_N_mdeg_hr_srHz / 60 * 1e-3 * sqrt(2 * pi);                     % [deg/srhr]

    % Extract errors from fit (Walker notation) (needs to be double checked, I did this in a rush)
    psd_N_err = (10^(psd_arw_fit.b + psd_arw_fit_conf) - 10^(psd_arw_fit.b)) / sqrt(2) / 60 * 1e-3 * sqrt(2 * pi);
    
    % Save quantities
    computed_quantities.psd_N_deg_srhr = psd_N_deg_srhr;
    computed_quantities.psd_N_deg_srhr_err = psd_N_err;
end


% RBI fit
psd_rbi_bool = all([f > psd_rbi_lim_s; f < psd_rbi_lim_e]);

if sum(psd_rbi_bool) > 1
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
    psd_rbi_fit_plot = 10.^feval(psd_rbi_fit, f_rbi_log) / 60;      % [mdeg/srhr]
    
    % Extract physical quantities from fit (IEEE notation)
    psd_B_mdeg_hr = 10^psd_rbi_fit.b * sqrt(2 * pi) / sqrt(2);      % [mdeg/hr]

    % Extract physical quantities from fit (Walker's notation)
    psd_B_deg_hr = psd_B_mdeg_hr * 1e-3;                                                % [deg/hr]

    % Extract errors from fit (Walker notation) (needs to be double checked, I did this in a rush)
    psd_B_err = (10^(psd_rbi_fit.b + psd_rbi_fit_conf) - 10^(psd_rbi_fit.b)) * sqrt(2 * pi) / sqrt(2) * 1e-3;
    
    % Save quantities
    computed_quantities.psd_B_deg_hr = psd_B_deg_hr;
    computed_quantities.psd_B_deg_hr_err = psd_B_err;
end


% % Generating linear fit data for plot
% psd_awn_fit_plot = 10.^feval(psd_awn_fit, f_awn_log) / 60;      % [mdeg/srhr]
% psd_arw_fit_plot = 10.^feval(psd_arw_fit, f_arw_log) / 60;      % [mdeg/srhr]
% psd_rbi_fit_plot = 10.^feval(psd_rbi_fit, f_rbi_log) / 60;      % [mdeg/srhr]
% 
% 
% % Extract physical quantities from fit (IEEE notation)
% psd_PHI_mdeg_hr_Hz3_2 = 10^psd_awn_fit.b / 2 / pi / sqrt(2);    % [mdeg/hr/Hz^(3/2)]
% psd_N_mdeg_hr_srHz = 10^psd_arw_fit.b / sqrt(2);                % [mdeg/hr/srHz]
% psd_B_mdeg_hr = 10^psd_rbi_fit.b * sqrt(2 * pi) / sqrt(2);      % [mdeg/hr]
% 
% % Extract physical quantities from fit (Walker's notation)
% psd_PHI_deg_srHz = psd_PHI_mdeg_hr_Hz3_2 / 3600 * 1e-3 * (2 * pi) * sqrt(2 * pi);   % [deg/srHz] (in Walker, the units of this term is deg)
% psd_N_deg_srhr = psd_N_mdeg_hr_srHz / 60 * 1e-3 * sqrt(2 * pi);                     % [deg/srhr]
% psd_B_deg_hr = psd_B_mdeg_hr * 1e-3;                                                % [deg/hr]
% 
% % Extract errors from fit (Walker notation) (needs to be double checked, I did this in a rush)
% psd_PHI_err = (10^(psd_awn_fit.b + psd_awn_fit_conf) - 10^(psd_awn_fit.b)) / 2 / pi / sqrt(2) / 3600 * 1e-3 * (2 * pi) * sqrt(2 * pi);
% psd_N_err = (10^(psd_arw_fit.b + psd_arw_fit_conf) - 10^(psd_arw_fit.b)) / sqrt(2) / 60 * 1e-3 * sqrt(2 * pi);
% psd_B_err = (10^(psd_rbi_fit.b + psd_rbi_fit_conf) - 10^(psd_rbi_fit.b)) * sqrt(2 * pi) / sqrt(2) * 1e-3;

%% Fit Allan
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

% work with log-log adev
avar_deg2_sec2 = avar_clean / (3600^2);             % [deg^2/s^2]
adev_log = log10(sqrt(avar_deg2_sec2));             % [deg/s]
tau2_log = log10(tau1);

% Generating a Piecewise model
allen_dev_model_piecewise_parts = ones(3, length(tau1));


% Data division and Fit
% AWN fit
awn_bool = tau2_log < awn_lim;

if sum(awn_bool) > 1
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
    
    % Extract physical quantities from fit (IEEE notation)
    allan_PHI_deg_srHz = 10^awn_allan_fit.b / sqrt(3) / sqrt(bw);       % [deg/srHz]

    % Extract physical quantities for piecewise
    allan_awn_deg = 10^awn_allan_fit.b;                                 % [deg]         PHI*sr(fH)

    % Extract physical quantities from fit (Walker notation)
    allan_PHI_deg_srHz = 10^awn_allan_fit.b / sqrt(3) / sqrt(bw);           % [deg/srHz] (in Walker, the units of this term is deg)

    % Extract errors from fit (Walker notation) (needs to be double checked, I did this in a rush)
    allan_PHI_err = (10^(awn_allan_fit.b + awn_allan_fit_conf) - 10^(awn_allan_fit.b)) / sqrt(3) / sqrt(bw);
    
    % Generating a Piecewise model
    allen_dev_model_piecewise_parts(1, :) = allen_dev_model_piecewise_parts(1, :) .* allan_awn_deg ./ tau1;                     % [deg/s]
    
    % Save quantities
    computed_quantities.allan_PHI_deg_srHz = allan_PHI_deg_srHz;
    computed_quantities.allan_PHI_deg_srHz_err = allan_PHI_err;
end

% ARW fit
arw_bool = all([tau2_log > arw_lim_s ; tau2_log < arw_lim_e]);

if sum(arw_bool) > 1
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
    
    % Extract physical quantities from fit (IEEE notation)
    allan_N_deg_s_srHz = 10^arw_allan_fit.b;                            % [deg/s/srHz] = [deg/srs]

    % Extract physical quantities for piecewise
    allan_arw_deg_s_srHz = 10^arw_allan_fit.b;                          % [deg/s/srHz]

    % Extract physical quantities from fit (Walker notation)
    allan_N_deg_srhr = 10^arw_allan_fit.b * 60;                             % [deg/srhr]

    % Extract errors from fit (Walker notation) (needs to be double checked, I did this in a rush)
    allan_N_err = (10^(arw_allan_fit.b + arw_allan_fit_conf) - 10^(arw_allan_fit.b)) * 60;
    
    % Generating a Piecewise model
    allen_dev_model_piecewise_parts(2, :) = allen_dev_model_piecewise_parts(2, :) .* allan_arw_deg_s_srHz ./ (tau1.^(0.5));     % [deg/s]
    
    % Save quantities
    computed_quantities.allan_N_deg_srhr = allan_N_deg_srhr;
    computed_quantities.allan_N_deg_srhr_err = allan_N_err;
end

% RBI fit
rbi_bool = tau2_log > rbi_lim;

if sum(rbi_bool) > 1
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
    allan_B_deg_s = 10^rbi_allan_fit.b * sqrt(pi / (2 * log(2)));       % [deg/s]

    % Extract physical quantities for piecewise
    allan_rbi_deg_s = 10^rbi_allan_fit.b;                               % [deg/s]

    % Extract physical quantities from fit (Walker notation)
    allan_B_deg_hr = 10^rbi_allan_fit.b * sqrt(pi / (2 * log(2))) * 3600;   % [deg/hr]

    % Extract errors from fit (Walker notation) (needs to be double checked, I did this in a rush)
    allan_B_err = (10^(rbi_allan_fit.b + rbi_allan_fit_conf) - 10^(rbi_allan_fit.b)) * sqrt(pi / (2 * log(2))) * 3600;
    
    % Generating a Piecewise model
    allen_dev_model_piecewise_parts(3, :) = allen_dev_model_piecewise_parts(3, :) .* allan_rbi_deg_s;                           % [deg/s]
    
    % Save quantities
    computed_quantities.allan_B_deg_hr = allan_B_deg_hr;
    computed_quantities.allan_B_deg_hr_err = allan_B_err;
end

% % Extract physical quantities from fit (IEEE notation)
% allan_PHI_deg_srHz = 10^awn_allan_fit.b / sqrt(3) / sqrt(bw);       % [deg/srHz]
% allan_N_deg_s_srHz = 10^arw_allan_fit.b;                            % [deg/s/srHz] = [deg/srs]
% allan_B_deg_s = 10^rbi_allan_fit.b * sqrt(pi / (2 * log(2)));       % [deg/s]
% 
% % Extract physical quantities for piecewise
% allan_awn_deg = 10^awn_allan_fit.b;                                 % [deg]         PHI*sr(fH)
% allan_arw_deg_s_srHz = 10^arw_allan_fit.b;                          % [deg/s/srHz]
% allan_rbi_deg_s = 10^rbi_allan_fit.b;                               % [deg/s]
% 
% % Extract physical quantities from fit (Walker notation)
% allan_PHI_deg_srHz = 10^awn_allan_fit.b / sqrt(3) / sqrt(bw);           % [deg/srHz] (in Walker, the units of this term is deg)
% allan_N_deg_srhr = 10^arw_allan_fit.b * 60;                             % [deg/srhr]
% allan_B_deg_hr = 10^rbi_allan_fit.b * sqrt(pi / (2 * log(2))) * 3600;   % [deg/hr]
% 
% % Extract errors from fit (Walker notation) (needs to be double checked, I did this in a rush)
% allan_PHI_err = (10^(awn_allan_fit.b + awn_allan_fit_conf) - 10^(awn_allan_fit.b)) / sqrt(3) / sqrt(bw);
% allan_N_err = (10^(arw_allan_fit.b + arw_allan_fit_conf) - 10^(arw_allan_fit.b)) * 60;
% allan_B_err = (10^(rbi_allan_fit.b + rbi_allan_fit_conf) - 10^(rbi_allan_fit.b)) * sqrt(pi / (2 * log(2))) * 3600;
% 
% % Generating a Piecewise model
% allen_dev_model_piecewise_parts = ones(3, length(tau1));
% allen_dev_model_piecewise_parts(1, :) = allen_dev_model_piecewise_parts(1, :) .* allan_awn_deg ./ tau1;                     % [deg/s]
% allen_dev_model_piecewise_parts(2, :) = allen_dev_model_piecewise_parts(2, :) .* allan_arw_deg_s_srHz ./ (tau1.^(0.5));     % [deg/s]
% allen_dev_model_piecewise_parts(3, :) = allen_dev_model_piecewise_parts(3, :) .* allan_rbi_deg_s;                           % [deg/s]


% Getting the piecewise of allan model fit
adev_deg_h_piecewise = max(allen_dev_model_piecewise_parts, [], 1) * 3600;      % [deg/h]

%% Prepare NGC curve
% calculate inertial rotation Allan Var % [sigy]=(deg/hr), [tau]=s
try 
    [sigy, tau] = (ADEV(Omega_deg_hr(2:end), fs)); 
catch; tau = nan; sigy = nan;
end

% NGC specs (Spin-Exchange Pumped NMR Gyros - arXiv:1604.03982)
ngc_awn = 0.012;    % [deg]
ngc_arw = 0.005;    % [deg/srhr]
ngc_bias = 0.02;    % [deg/hr]

% Generate the NGC PSD and ADEV lines

% Normalizing to world rotation pure units [Hz^2/Hz]
ngc_h2 = (2 * pi * ngc_awn / 360 / sqrt(3 * fH_ngc))^2;                 % [1/Hz] 
ngc_h0 = (sqrt(2) * ngc_arw * sqrt(2 * pi) / 360 / sqrt(3600))^2;       % [Hz]
ngc_hm1 = (ngc_bias / sqrt(2 * log(2)) / 360 / 3600)^2;                 % [Hz^2]

% Generating NGC spec PSD line in sqrt([Hz^2/Hz]) units
ngc_psd = sqrt(ngc_h2 * f.^2 + ngc_h0 * f.^0  + ngc_hm1 * f.^(-1));     % [Hz^2/Hz]
ngc_psd_mdeg_srhr = ngc_psd * 1e3 / sqrt(2 * pi) * 360 * sqrt(3600);

% Allan Deviation of NGC spec
ADAV_NGC = sqrt((ngc_awn ./ (tau / 3600)).^2 + (ngc_arw ./ sqrt(tau / 3600)).^2 + ngc_bias^2);
% ADAV_gen = sqrt((PHI_deg_srhr_Hz ./ (tau / 3600)).^2 + (N_deg_srhr ./ sqrt(tau / 3600)).^2 + B_deg_rthr_rtHz^2); % / sqrt(360)

%% Plotting the Data
% Setting up the legend text
awn_text = ['AQN/AWN (\Phi:f_H=' num2str(fit_bw,2) ') [^\circ/\surd{Hz}]', ' | ' 'Allan: ', num2str(allan_PHI_deg_srHz), ' | ', 'Psd: ', num2str(psd_PHI_deg_srHz), ' | '];
arw_text = ['ARW(N) [^\circ/\surd{hr}]', '  |           ' 'Allan: ', num2str(allan_N_deg_srhr), ' | ', 'Psd: ', num2str(psd_N_deg_srhr), ' | '];
rbi_text = ['RBI(B)   [^\circ/hr]  ', '  |           ' 'Allan: ', num2str(allan_B_deg_hr), ' | ', 'Psd: ', num2str(psd_B_deg_hr), ' | '];


figure1 = gcf; 
set(figure1, 'Position',[1 41 1920 963])

% Plot measured data
axes1 = axes('Parent', figure1, 'Position', [0.054 0.6843 0.4073 0.2741]);
hold(axes1,'on');
yyaxis right
plot(t,(dBz-mean(dBz))*1e6*1e3);
% plot(t,(dBz-mean(dBz))*1e6*1e3,'color',[0 0 0.8]);
ylabel('\DeltaB_z[nG]','fontsize',16);
yyaxis left
plot(t, Omega_deg_hr - mean(Omega_deg_hr));
% plot(t, Omega_deg_hr - mean(Omega_deg_hr),'color',[0.6 0.2 0]);
legend('\Omega [^\circdeg/hr]','\DeltaB_z[nG] ');
grid(axes1, 'on');
set(axes1, 'FontSize',16);
xlabel('Time (sec)','fontsize',16);
ylabel('\Omega [^\circdeg/hr]','fontsize',16);
hold off

% Plot PSD
axes2 = axes('Parent', figure1, 'Position', [0.054 0.3531 0.6703 0.2368]);
hold(axes2,'on');
xlabel('Frequency (Hz)', 'fontsize', 16); 
ylabel('\delta\omega (mdeg/\surdhr)', 'fontsize', 16); 
loglog(f, N_mdeg_rthr, 'color', [0 0.4470 0.7410], 'DisplayName', '\Omega');
plot(f, ngc_psd_mdeg_srhr, ':k', 'linew', 3, 'DisplayName', 'NGC');

% Plot PSD fit
if sum(psd_awn_bool) > 1
    plot(f_awn_linear, psd_awn_fit_plot, ':g', 'linew', 2, 'DisplayName', awn_text)
end

if sum(psd_arw_bool) > 1
    plot(f_arw_linear, psd_arw_fit_plot, ':r', 'linew', 2, 'DisplayName', arw_text)
end

if sum(psd_rbi_bool) > 1
    plot(f_rbi_linear, psd_rbi_fit_plot, ':m', 'linew', 2, 'DisplayName', rbi_text)
end

box(axes2,'on');
grid(axes2,'on');
axis(axes2,'tight');
axes2.YAxis.FontSize = 16;
axes2.XAxis.FontSize = 16;
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
% lgnd1=get(gca,'Children');
legend2 = legend;%([lgnd1(5),lgnd1(4)]);% display only Omega and NGC
set(legend2, 'Position', [0.2296 0.0765 0.3057 0.1796]);
set(legend2,'fontsize', 16);
% set(legend2,'Location', 'Best');

% Plot DBz

% axes2 = axes('Parent', figure1, 'Position', [0.054 0.3531 0.6703 0.2368]);
% hold(axes2,'on');
% hold on
% yyaxis right
% xlabel('Frequency (Hz)', 'fontsize', 16); 
% ylabel('\delta\omega (mdeg/\surdhr)', 'fontsize', 16); 
% loglog(f, N(1,:), 'color', [0.9 0.1 0.17410], 'DisplayName', 'dBz');
% 


% Plot Allan Deviation
axes3 = axes('Parent', figure1, 'Position', [0.5724 0.6854 0.3927 0.2772]);
loglog(tau1, sqrt(avar_clean), 'linew', 2);
hold(axes3,'on');
plot(tau, ADAV_NGC, '--k', 'linew', 2)

if sum(awn_bool) > 1
    plot(tau1(awn_bool), allen_dev_model_piecewise_parts(1, awn_bool) * 3600, ':g', 'linew', 2)
end
if sum(arw_bool) > 1
    plot(tau1(arw_bool), allen_dev_model_piecewise_parts(2, arw_bool) * 3600, ':r', 'linew', 2)
end
if sum(rbi_bool) > 1
    plot(tau1(rbi_bool), allen_dev_model_piecewise_parts(3, rbi_bool) * 3600, ':m', 'linew', 2)
end

% % sanity check for ADEV
% aaaaa = ngc_awn * ones(size(tau1)) * 3600;
% plot(tau1(awn_bool), aaaaa(awn_bool)./tau1(awn_bool), ':r', 'linew', 2, 'Parent', axes3)

grid on;
xlabel('\tau (s)', 'fontsize',16);
ylabel('\sigma_{\ity} (deg/hr)', 'fontsize',16); 
% title('Allan deviation');
axes3.YAxis.FontSize = 16;
axes3.XAxis.FontSize = 16;


% mkNicePlt;

warning on
warning('Need to divide ARW from PSD by 2.4 to get the same as in simulation.')
end

