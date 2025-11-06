% Data for Yogev
% The following code extracts the measured signals and scales them relative
% to input, namely relative to the output of TIA.
load('experiment_raw_data.mat');
analyze_ESR = 1;
analyze_NMR_129 = 0;
analyze_NMR_131 = 0;
analyze_TIA = 0;

% Plot RAW data
% figure(3452), clf
% plot(meas.cDAQ_t,meas.cDAQ_v);
% legend(meas.lgd_cell);
% xlabel('Time [s]');
% ylabel('Amplitude [V]');


% Fs = prm.cDAQ_Fs; %sampling frequency Sa/s
sensitivity_ESR = 0.5; % the sensitivity of ESR LIA
sensitivity_NMR_129 = prm.WP.LIAconfig.SR129.Sensitivity; % the sensitivity of NMR LIA 129
sensitivity_NMR_131 = prm.WP.LIAconfig.SR131.Sensitivity; % the sensitivity of NMR LIA 129
ga = (-7.44e7  )/(2*pi)/(1e4);              % [Hz/G] Xe129 gyromagnetic ratio
gb = (+2.2056e7)/(2*pi)/(1e4);              % [Hz/G] Xe131 gyromagnetic ratio
x = ga / gb;

% PID calibarations from voltage to magnetic field and rotation rate
dDeg_per_hr_dV = 360 * 3600*prm.WP.modulationFactors.df131/ 5 / (1 - 1 / x); %deg/hr / Volt
dBz_dV = 320e-6; % Gauss/ Volt

% gain of LIA NMR
%  The following gains gives the following signal: output = input*gain. In
%  order to translate a measured output signal to relative to input we need
%  to divide by the gain.
voltage_divider = 7.89;
ESR_theta_atten_factor = 10;
NMR_131_X_and_R_atten_factor = 5.3;
LIA_ESR_gain = 10/sensitivity_ESR/voltage_divider;
LIA_NMR_129_gain = LIA_ESR_gain/voltage_divider*10/sensitivity_NMR_129;
LIA_NMR_131_gain = LIA_ESR_gain/voltage_divider*10/sensitivity_NMR_131;


% Definition of the variables
if analyze_ESR
%     ESR_LIA_R_volt = meas.cDAQ_v(:,1)/LIA_ESR_gain;
    ESR_LIA_X_volt = meas.cDAQ_v(:,1)/LIA_ESR_gain;
    ESR_LIA_Y_volt = meas.cDAQ_v(:,2)/LIA_ESR_gain;
%     ESR_LIA_Theta_volt = meas.cDAQ_v(:,4);
%     ESR_LIA_Theta_deg = meas.cDAQ_v(:,4)/10*180*ESR_theta_atten_factor;
end
if analyze_NMR_129
%     LIA_129_Y_volt = meas.cDAQ_v(:,8)/LIA_NMR_129_gain; % Relative to TIA
%     LIA_129_X_volt = meas.cDAQ_v(:,7)/LIA_NMR_129_gain; % Relative to TIA
    LIA_129_R_volt = meas.cDAQ_v(:,3)/LIA_NMR_129_gain; % Relative to TIA
    LIA_129_PID_volt = meas.cDAQ_v(:,5);
    LIA_129_Theta_volt = meas.cDAQ_v(:,2);
    LIA_129_PID_Gauss = LIA_129_PID_volt*dBz_dV;
    LIA_129_Theta_deg = LIA_129_Theta_volt* 180/10;
end
if analyze_NMR_131
%     LIA_131_Y_volt = meas.cDAQ_v(:,10)/LIA_NMR_131_gain; % Relative to TIA
%     LIA_131_X_volt = meas.cDAQ_v(:,9)/LIA_NMR_131_gain*NMR_131_X_and_R_atten_factor; % Relative to TIA
    LIA_131_R_volt = meas.cDAQ_v(:,4)/LIA_NMR_131_gain*NMR_131_X_and_R_atten_factor; % Relative to TIA
    LIA_131_PID_volt = meas.cDAQ_v(:,1);
    LIA_131_Theta_volt = meas.cDAQ_v(:,6);
    LIA_131_PID_deg_per_hr = LIA_131_PID_volt*dDeg_per_hr_dV;
    LIA_131_Theta_deg = LIA_131_Theta_volt* 180/10;
end
if analyze_TIA
   meas.time_sec = meas.scope_experiment.time;
   meas.TIA_volt = meas.scope_experiment.cha1(1,:);
end

%% Plot ESR
if analyze_ESR
    figure(3452), clf
    tiledlayout('flow')
%     % R
%     nexttile
%     plot(meas.cDAQ_t,ESR_LIA_R_volt);
%     xlabel('Time [s]');
%     ylabel(['ESR LIA R (Volt RMS)']);
%     grid on
    
    % X
    nexttile
    plot(meas.cDAQ_t,ESR_LIA_X_volt);
    xlabel('Time [s]');
    ylabel(['ESR LIA X (Volt RMS)']);
    grid on
    
    % Y
    nexttile
    plot(meas.cDAQ_t,ESR_LIA_Y_volt);
    xlabel('Time [s]');
    ylabel([ 'ESR LIA Y (Volt RMS)']);
    grid on
    
    % Theta
%     nexttile
%     yyaxis left
%     plot(meas.cDAQ_t,ESR_LIA_Theta_volt);
%     xlabel('Time [s]');
%     ylabel(['ESR LIA Theta (Volt)']);
%     grid on
%     yyaxis right
%     plot(meas.cDAQ_t,ESR_LIA_Theta_deg);
%     ylabel([ 'ESR LIA Theta (Degrees)']);
end

%% Plot NMR 129
if analyze_NMR_129
    figure(3452), clf
    tiledlayout('flow')
%     % Y
%     nexttile
%     plot(meas.cDAQ_t,LIA_129_Y_volt);
%     xlabel('Time [s]');
%     ylabel(['LIA 129 Y (Volt RMS RMS)']);
%     grid on
    
    % X
%     nexttile
%     plot(meas.cDAQ_t,LIA_129_X_volt);
%     xlabel('Time [s]');
%     ylabel(['LIA 129 X (Volt RMS RMS)']);
%     grid on
    
    % R
    nexttile
    plot(meas.cDAQ_t,LIA_129_R_volt);
    xlabel('Time [s]');
    ylabel(['LIA 129 R (Volt RMS RMS)']);
    grid on
    
    % Theta
    nexttile
    yyaxis left
    plot(meas.cDAQ_t,LIA_129_Theta_volt);
    xlabel('Time [s]');
    ylabel([ 'LIA 129 Theta (Volt)']);
    grid on
    yyaxis right
    plot(meas.cDAQ_t,LIA_129_Theta_deg);
    ylabel(['LIA 129 Theta (Degrees)']);
    
    % PID out
    nexttile
    yyaxis left
    plot(meas.cDAQ_t,LIA_129_PID_volt);
    xlabel('Time [s]');
    ylabel(['PID 129 out (Volt)']);
    grid on
    yyaxis right
    plot(meas.cDAQ_t,LIA_129_PID_Gauss);
    ylabel(['PID 129 out (Gauss)']);
    
end

%% Plot NMR 131
if analyze_NMR_129
    figure(3453), clf
    tiledlayout('flow')
%     % Y
%     nexttile
%     plot(meas.cDAQ_t,LIA_131_Y_volt);
%     xlabel('Time [s]');
%     ylabel(['LIA 131 Y (Volt RMS RMS)']);
%     grid on
    
%     % X
%     nexttile
%     plot(meas.cDAQ_t,LIA_131_X_volt);
%     xlabel('Time [s]');
%     ylabel(['LIA 131 X (Volt RMS RMS)']);
%     grid on
    
    % R
    nexttile
    plot(meas.cDAQ_t,LIA_131_R_volt);
    xlabel('Time [s]');
    ylabel([ 'LIA 131 R (Volt RMS RMS)']);
    grid on
    
    % Theta
    nexttile
    yyaxis left
    plot(meas.cDAQ_t,LIA_131_Theta_volt);
    xlabel('Time [s]');
    ylabel(['LIA 131 Theta (Volt)']);
    grid on
    yyaxis right
    plot(meas.cDAQ_t,LIA_131_Theta_deg);
    ylabel(['LIA 131 Theta (Degrees)']);
    
    % PID out
    nexttile
    yyaxis left
    plot(meas.cDAQ_t,LIA_131_PID_volt);
    xlabel('Time [s]');
    ylabel(['PID 131 out (Volt)']);
    grid on
    yyaxis right
    plot(meas.cDAQ_t,LIA_131_PID_deg_per_hr);
    ylabel(['PID 131 out (deg/hr)']);
    
end

%% Generate PSD and allan for 131
if analyze_NMR_131
    addpath('src')
    addpath('utils')
    addpath('instruments')
    
    meas.time_s = meas.cDAQ_t;
    meas.Det131_PIDb_V = smooth(LIA_131_PID_volt, prm.data_smooth_param);LIA_131_PID_volt;%;%;
    meas.Bz_V = smooth(LIA_129_PID_volt, prm.data_smooth_param);%meas.Det131_PIDb_V;%smooth(LIA_131_PID_volt, prm.data_smooth_param);LIA_131_PID_volt;%;%; % so the code does not crash
    prm.df131 = prm.WP.modulationFactors.df131;
    
    [f,N,T,S,tau, noise_mdl_prm] = get_noise_spectrum_with_fit(meas, prm, prm.estimated_bw);
    meas.noise_mdl_prm = noise_mdl_prm;
    %     plot_single_experiment(meas, prm)
end

%% Plot TIA
if analyze_TIA
    figure(235434), clf
plot(meas.time_sec, meas.TIA_volt,'-');
xlabel('Time [s]');
ylabel('TIA out [V]');
end
