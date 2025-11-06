% close all

%% How to work with the code
% The flags are stored in a flg struct. The parameters that are saved after
% the run are saved in prm structure. All the measurements are stored in a
% meas structure. In the end, only those structures are saved. All other
% variables are not saved and should be derived given the saved structures.
%

%% TODO

% add amplitude reading for Rb FID (Kostya), zero time shift. 
% fix audiodata
% write instruments map (structure)
% write Xenon FID experiment
% connect 1/2 lambda to probe
% write a probe intensity PSD experiment
% add working point parameters to logger
% compute the error in T2 in the "FID rubidium" script
% check that fft1 does not use dc (we removed it)
warning('check that fft1 does not use dc (and we removed it from code)');

%% Reset ?
% reset_all = 1; % uncomment if you need to force reset
addpath('src')
addpath('utils')
addpath('instruments')
substitute_if_does_not_exist('reset_all',1); % clears all variables and will reconnect to instruments
clear meas flg prm % clear vars from prev measuremnts
if reset_all
    clear all;
    reset_all = 1;
    instrreset;
    % add to path list

end

%% Experiment description
exp_description = ['NMR LIA R theta + PID 129+131 2kSps 60min']; 

%% Flags

% experiments
flg.measure_world_rotation =                  0;
flg.FID_Rb =                                  0;                             % 
flg.Rb_absorption_vs_Laser_wavelength =       0;         % 0 - no measure, 1 - measure vs Pump laser, 2 - measure vs Prob laser
flg.FID_Rb_vs_Laser_wavelength =              0;         % 0 - no measure, 1 - measure vs Pump laser, 2 - measure vs Prob laser, 3 - measure vs Pump and Probe laser
flg.FID_Rb_vs_Pump_And_Probe_detuning =       0;  % 0 - no measure, 1- regular measure, 2- Noam's measure
flg.photon_shot_noise =                       0;
flg.FID_Xenon =                               0;
flg.Compute_PSD =                             0;
flg.get_data_from_scope =                     0;
flg.get_data_from_cDaq =                      1; %(after measurement, data is available in the workspace in the "meas" struct)
flg.lia_noise_floor =                         0;
flg.load_past_experiment =                    0;

% options
% load working point
flg.load_wp = 0;
flg.load_laser_params = 0;                  % use 0 to leave laser detunings as is
flg.plot_data = 1;
flg.close_mid_plots = 0;                    % closing all plots which are not generated at end of experiment

% save data
flg.save_data = 1;                       

%% Input Parameters
% get useful constants
[unt, P, c] = units_consts;

% world roatation measurement
flg.get_noise_spectrum = 1;                 % only enabled with world_rot_meas
prm.exp_time_min = 60;
prm.points_to_plot = 1000;
prm.data_smooth_param = 100;                % [100]
prm.estimated_bw = 3.0;                     % [Hz]

% FID Rb / Rb absorption vs pump wavelength measurement
prm.pump_span_volts = [-0.02, 0.02];
prm.pump_span_points = 21;

% FID Rb / Rb absorption vs probe wavelength measurement
prm.probe_span_volts = [0, 0];
prm.probe_span_points = 10;

% Photon Shot Noise
prm.psn_measurements = 16; % pick from [1, 4, 9, 16, 25, ...]
prm.psn_reps = 3;

% Get data from scope
prm.channels_scope = [1];
% gain = 10 / 20e-3;
% under_sample = 100000;
prm.gains_scope = [1 1 1 1];
% set scope for Xenon FID measurement
% Sampling rate depends on the number of active channels. The scope has
% maximal memory of 62500 samples (when all channels are active without high res). Sampling
% rate in this case is Fs = 62500/Tscale/10;
prm.Fs_scope = 12.5e6; % Sp/sec
prm.Tscale_scope = 62500/prm.Fs_scope/10; % sec/div


% Data Acquisition directly from NI cDAQ-9174                                
% (card 1 or 2, channel 0 to3 ) [(2,0)  (2,1) (1,0) (1,1) (2,2) (2,3) (1,2) (1,3)  (3, 0) (3, 1) (3, 2) (3, 3)]
% prm.cDAQ_channels =           [  1      2     3     4     5     6     7     8      9      10     11     12 ];
prm.cDAQ_channels = [2 3 5 6 7 8]; % A vector of the channels we want to record;
                       %Use cDAQ.getChannels() to get a list of available channels
prm.cDAQ_meas_time_min = 60;
prm.cDAQ_Fs = 2048; %samples/sec (instrument default is 2048; max is 51200)

% Get PSD of signals
variable_time = "meas.scope_experiment.time";
variable_names = ["meas.scope_experiment.cha2"];

prm.DAQ_ATT_129 = 2.22;   % DAQ Xe129 phase attenuation factor
prm.DAQ_ATT_131 = 2.25;   % DAQ Xe131 phase attenuation factor
prm.PAUSE_TIME_SEC = 1;   % pause after measurement

% Paths
save_path = 'D:\NMRG\NMRG_code_and_data_history';

%% Initialize

if exist('mult_flg_is_running','var') && mult_flg_is_running
    % Multiple measurements:
    % Automatically change the experiment description
    exp_description = [exp_description mult_cur_exp_description];
end

if ~flg.save_data
    warning on
    warning('The data will not be saved.')
end

gongData = load('gong.mat');
gongObj = audioplayer(gongData.y,gongData.Fs);

% get main script path
main_path = pwd;


meas.DateTime_var = datetime;
meas.DateNumber = datenum(meas.DateTime_var);

% communication with instruments
if reset_all
    connect_to_instruments_Itay
    get_daq
end

%Initial Logger
logger = ExperimentLogger(fullfile(main_path, 'log'));
if 1
    logger.log('#########################################################')
    logger.log('Experiment Description:')
    logger.log(exp_description)
    logger.log('#########################################################')
    disp('#########################################################')
    disp('Experiment Description:')
    disp(exp_description)
    disp('#########################################################')
end

%% Read / Load working point
 
if flg.load_wp
    load_working_point;
else
    read_current_working_point;
end

%% Read Scope
if flg.get_data_from_scope
   get_data_from_scope 
end

%% Read cDAQ
if flg.get_data_from_cDaq
   get_data_from_cDAQ 
end

%% Measure FID rubidium
if flg.FID_Rb
    % ========================================================================
    %   This script is measuring Alkali FID and computing the T2 and resonance 
    %   frequncies from the data. The code assumes Bz is set, Bx and By are 
    %   zeroed, BPD is balanced by a lambda plate. 
    %   The signal from BPD is connected ch1 scope1. 
    % ========================================================================
    try
        Alkali_FID_measurement;
    catch e
        disp(['a problem occured in Alkali_FID_measurement. The message was: ' e.message])
        logger.log(['a problem occured in Alkali_FID_measurement. The message was: ' e.message])
    end
    
end

%% Measure FID rubidium vs Laser wavelength
if flg.FID_Rb_vs_Laser_wavelength
    % ========================================================================
    %   This script is used for Alkali FID measurement vs laser detunings.
    %   In order to choose between the pump / probe lasers, set the experiment
    %   flag flg.FID_Rb_vs_Laser_wavelength to 1 / 2 respectively.
    %
    %   This script uses the Alkali_FID_measurement script, and based on the 
    %   assumptions written in the headline of the that script.
    % ========================================================================
    Alkali_FID_vs_laser_detuning
    
end

%% Measure rubidium absorption vs Laser wavelength
if flg.Rb_absorption_vs_Laser_wavelength
    % ========================================================================
    %   This script is similar to Alkali FID measurement vs laser detunings.
    %   In order to choose between the pump / probe lasers, set the experiment
    %   flag flg.Rb_absorption_vs_Laser_wavelength to 1 / 2 respectively.
    %
    %   Note: We don't Disable channels 2 of AG1 and AG2! the user should
    %   set them to the desired values.
    %   To measure the absorption connect the PD diff, -, and + outputs to
    %   channels 1, 3, and 4 (respectively) of scope 1.
    %   Note: the PD diff should already be connected to ch1 of scope 1.
    % ========================================================================
    Alkali_absorption_vs_laser_detuning
    
end

%% Measure FID rubidium vs both Lasers detunings

if flg.FID_Rb_vs_Pump_And_Probe_detuning==1
    % ========================================================================
    %   This script is used for Alkali FID measurement vs probe laser and pump 
    %   laser detunings simultaneously.
    %   The script sets a probe detuning value and then sweeps the pump
    %   detuning for that value.
    %
    %   This script uses the next scripts:
    %       Alkali_FID_measurement.m 
    %       Alkali_FID_vs_laser_detuning.m
    %   and based on the assumptions written in the headlines of those scripts.
    % ========================================================================
    Alkali_FID_vs_probe_and_pump_detuning;
end
if flg.FID_Rb_vs_Pump_And_Probe_detuning==2
    FID_Good_Index=0;
    FIDindex=1;
    
    Alkali_FID_vs_probe_and_pump_detuning_noam;
end

%% Measure World Rotation
if flg.measure_world_rotation
    % ========================================================================
    %   This script implementing a World Rotation Measurement and computing the 
    %   power spectral density (PSD) from the measured data. It assumes that 
    %   the system is in an optimal working point.
    % ========================================================================
    measure_world_rotation;
    
end

%% Measure Photon Shot Noise
if flg.photon_shot_noise
    % =====================================================================
    %   This script is used for Photon Shot Noise measurement.
    %
    %
    % =====================================================================
    
    photon_shot_noise_measurement;
end

%% Measure Xenon FID
if flg.FID_Xenon
    % =====================================================================
    %   This script is for measuring Xenon FID, computing the T2s and 
    %   resonance frequncies from the data. The code assumes Bz, Bx and By 
    %   are set, NMRs are working and the PIDs are disabled. 
    %   The signal is measured from the scope's channel 2 which is
    %   assumed to be the signal after the ESR LIA.
    % =====================================================================
    try
        Xenon_FID;
    catch e
        disp(['a problem occured in Xenon FID measurement. The message was: ' e.message])
        logger.log(['a problem occured in Xenon FID measurement. The message was: ' e.message])
    end
end

%% Compute the Power Spectral Density of data
    % =====================================================================
    %   Here we compute the Power Spectral Density of the signals given in
    %   variable_names with fs computed from variable_time.
    % =====================================================================

if flg.Compute_PSD
    data2psd = [];
    names2psd = [];
    
    % get the signals
    for i = 1:length(variable_names)
        if ~exist(variable_names(i), 'var')
            meas_variable = eval(strcat(variable_names(i), ';'));
            data2psd = [data2psd meas_variable'];
            names2psd = [names2psd variable_names(i)];                   
        end
    end
    
    % get the time vector
    time_variable = eval(strcat(variable_time, ';')); 
    fs2psd = 1 / mean(diff(time_variable));
    
    % compute the PSD
    get_noise_spectrum(data2psd, fs2psd, names2psd);
end

%% LIA noise floor measurement
if flg.lia_noise_floor
    % =====================================================================
    %   This script is for measuring the white noise floor of the LIA with
    %   50 ohm termination at the input
    % =====================================================================
    lia_noise_floor_measurement
end

%% Load previous experiment results and parameters

if flg.load_past_experiment
    % =====================================================================
    %   Loading results of an old experiment 
    % =====================================================================
    address = 'D:\NMRG\NMRG_code_and_data_history\2022';
    [fname, address] = uigetfile(address, 'Select an experiment to load to system');
    temp = load([address fname], 'prm', 'meas', 'flg');
    prm = temp.prm;
    meas = temp.meas;
    flg = temp.flg;
    
    % We do not want to resave the old results again
    flg.save_data = 0;
    
end

%% Save Data
logger.log('Done')

if flg.save_data
    logger.log('Saving everything')
    save_everything
end

%% Finish experiment notification
reset_all = 0;  % set reset to zero
playblocking(gongObj);