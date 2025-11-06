% close all

%% How to work with the code
% The flags are stored in a flg struct. The parameters that are saved after
% the run are saved in prm structure. All the measurements are stored in a
% meas structure. In the end, only those structures are saved. All other
% variables are not saved and should be derived given the saved structures.
%

%% TODO
% Script copied from minime. MAKE EVERYTHING WORK IN DOR B

% add windowing to PSD analysis (NIST book page 70)
% add data averaging for longer measuremnts (NIST book page 70)
% add amplitude reading for Rb FID (Kostya), zero time shift. 
% write instruments map (structure) so that we can use AWG interchangebly
% write a probe intensity PSD experiment
% add working point parameters to logger
% compute the error in T2 in the "FID rubidium" script, correct the problem
% with fitting to heating frequency istead of the second rb
% check that fft1 does not use dc (we removed it)
% add each measuremnt name to exp description
% warning('check that fft1 does not use dc (and we removed it from code)');
% connect Tabor heater to tcp ip
% Use numerr2str  to display results with error bars correctly!!!!!!!!!!!!!
% Fix runs to new standards of logger, loops and GUI

%% Reset ?
% reset_all = 1; % uncomment if you need to force reset
% [filepath,~,~] = fileparts(matlab.desktop.editor.getActiveFilename);
filepath = 'D:\Synced team folder\SynologyDrive\Dor B\DorB_operational';
cd(filepath); % Just in case
addpath('src')
addpath(genpath('utils'))
addpath('instruments')
substitute_if_does_not_exist('reset_all',1); % clears all variables and will reconnect to instruments
clear meas flg prm % clear vars from prev measuremnts
if reset_all
    clear
    reset_all = 1;
    instrreset;
    % add to path list
end

%% Experiment description
exp_description = ['']; % 

%% Loop data
substitute_if_does_not_exist('is_part_of_loop',0);

%% Flags

% experimentsse 
flg.measure_world_rotation =                  1; % OK 
flg.FID_Rb =                                  0; % OK
flg.FID_Xenon =                               0; % OK
flg.find_ESR_WP =                             0; % OK
flg.measure_mgnt_sens =                       0; % 

flg.Alkali_polarization =                     0; % OK need to install AOM and PD

flg.get_data_from_scope =                     0; % OK DorB        %Record data drom scope. See configurations below
flg.get_data_from_cDaq =                      0; % OK DorB       (after measurement, data is available in the workspace in the "meas" struct)

% flg.Rb_absorption_vs_Laser_wavelength =       0;         % 0 - no measure, 1 - measure vs Pump laser, 2 - measure vs Prob laser
% flg.FID_Rb_vs_Laser_wavelength =              0; % OK DorB    % 0 - no measure, 1 - measure vs Pump laser, 2 - measure vs Prob laser, 3 - measure vs Pump and Probe laser
% flg.FID_Rb_vs_Pump_And_Probe_detuning =       0;              % 0 - no measure, 1- regular measure, 2- Noam's measure
% flg.photon_shot_noise =                       0;
% flg.Compute_PSD =                             0;
% flg.lia_noise_floor =                         0;
%Flags marked with "OK DorB" were verified to work in DOR B
%Flags marked with "OK" are also suitable for new standards of logger, loops and GUI

% options
% load working point
flg.load_wp =                                 0;
flg.load_laser_params =                       0; % use 0 to leave laser detunings as is
flg.plot_data =                               0; % plots raw data from DAQ
flg.close_mid_plots =                         0; % closing all plots which are not generated at end of experiment
flg.save_all_open_figs =                      0; % when 0 saves only figures that are created using create_fig
flg.is_part_of_loop = is_part_of_loop;

% save data
flg.save_data = 1;                       

%% Input Parameters
% get useful constants
[unt, P, c] = units_consts;
add_experimental_calibration;

% daq definitions
prm.cdaq_pause_time_sec = 1;
prm.exp_time_min = 40; 
prm.save_time_min = 1; 
prm.ch_to_save= [1:12]; 
% prm.ch_to_save= [1 2 3 4 5 6 7 8]; 

prm.Fs = 2048; %[samples/sec] (instrument default is 2048; max is 51200)
prm.data_smooth_param = 1 + 15 * (prm.exp_time_min > 4);             
prm.max_full_variable_size_GB = 2;
flg.allow_variables_above_max = 1;
flg.keep_raw_data = (prm.exp_time_min <= 4) & prm.data_smooth_param > 1;

% world rotation measurement (see also daq definitions)
flg.verify_pid_locks = 0;
prm.locking_Bz_instead_of_drive = [0 0]; % for 129 and 131 respectively 
flg.plot_measured_signals = 1;
flg.analyze_rotation_signal = 1;          %1 for closed-loop analysis 
prm.estimated_bw = 8;                     % [Hz]
flg.plot_individual_xenons = 1;
prm.base_fig_num = 300;

% legacy - to remove
prm.american_lock = 0;                      % 1= PID on both NMR drives, 0= PID129 on magnetic field 
prm.points_to_plot = 1000;
flg.get_noise_spectrum = 1;                 % only enabled with world_rot_meas


% Alkali polarization measurement
prm.num_voltage_samples = 11;

% Xe FID options
% if flg.FID_Xenon
    flg.measure_magnetometer_responsivity = 1; % 0 if conneted to FPGA
    flg.use_DAQ_4_Xe_FID = 0; % 0 if taking from scope
%     prm.cDAQ_channels = [1 ]; % A vector of the channels we want to record;
%     %Use cDAQ.getChannels() to get a list of available channels
%     prm.cDAQ_meas_time_min = 20/60;
%     prm.cDAQ_Fs = 2*2048; %samples/sec (instrument default is 2048; max is 51200)
% end

% Photon Shot Noise
prm.psn_measurements = 16; % pick from [1, 4, 9, 16, 25, ...]
prm.psn_reps = 3;


% Get data from scope
% From which scope object to record? (input the var name as string)
prm.scope_obj = 'scope2'; 
% What channels?
prm.channels_scope = [1 2];
% What gains to use for each channel?
prm.gains_scope = [1 1 1 1];


% if flg.get_data_from_cDaq 

% Data Acquisition directly from NI cDAQ-9174                                
% (card 1 or 2, channel 0 to3 ) [ (1,0) (1,1) (1,2) (1,3) (2,0)  (2,1) (2,2) (2,3)  (3, 0) (3, 1) (3, 2) (3, 3)]
% prm.cDAQ_channels =           [  1      2     3     4     5     6     7     8      9      10     11     12 ];
prm.cDAQ_channels = 7; [2 5 6 7 8 9 10]; % A vector of the channels we want to record;
                       %Use cDAQ.getChannels() to get a list of available channels
prm.cDAQ_meas_time_min = 5;
prm.cDAQ_Fs =2048; %samples/sec (instrument default is 2048; max is 51200)
% end

% Get PSD of signals
variable_time = "meas.scope_experiment.time";
variable_names = ["meas.scope_experiment.cha2"];

prm.DAQ_ATT_129 = 2.22;   % DAQ Xe129 phase attenuation factor
prm.DAQ_ATT_131 = 2.25;   % DAQ Xe131 phase attenuation factor
prm.PAUSE_TIME_SEC = 1;   % pause after measurement

if flg.is_part_of_loop 
    if isvarname('full_loop_data')
        prm.loop_data = full_loop_data;
    else
        prm.loop_data = struct;
    end
end

% magnetomectic sensativity and LIA noise floor
prm.mgnt_exp_time_s = 10;
prm.mgnt_points_to_plot = 5000;
prm.mgnt_Fs = 2*2048;                               %[samples/sec] (instrument default is 2048; max is 51200)
prm.mgnt_meas_magsens = 1; %measure magnetometric responsivity 1, just noise floor 0
turn_on_Xe = 0; % this turns on Xe drives after the end of the measuremnt


% Paths
save_path = 'D:\Synced team folder\Data Archive\DorB\All Data\Currently synced data';

%% Initialize

% if exist('mult_flg_is_running','var') && mult_flg_is_running
    % Multiple measurements:
    % Automatically change the experiment description
%     exp_description = [exp_description mult_cur_exp_description];
% end

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
    connect_to_instruments_Dor_B
    get_daq
end

%Initial Logger
logger = ExperimentLogger(fullfile(main_path, 'log'),[],~flg.is_part_of_loop);
if 1
    logger.log('#########################################################')
    logger.log('Experiment Description:')
    logger.log(exp_description)
    if is_part_of_loop
        logger.log(loop_parse_loop_data(prm.loop_data));
    end

    logger.log('#########################################################')
end

%% Read / Load working point
 
if flg.load_wp
    load_working_point;
    flg.save_data = 0;
else
    save_current_working_point;
end


%% Run measurements according to flags
exp_flgs = {...
'FID_Rb',...
'Alkali_polarization',...
'FID_Xenon',...
'measure_world_rotation',...
'measure_mgnt_sens',...
'get_data_from_cDaq',...
'get_data_from_scope',...
'find_ESR_WP',...
% 'Rb_absorption_vs_Laser_wavelength',...
% 'FID_Rb_vs_Laser_wavelength',...
% 'FID_Rb_vs_Pump_And_Probe_detuning',...
% 'photon_shot_noise',...
% 'Compute_PSD',...
% 'lia_noise_floor',...
};
for ind_flg = 1:length(exp_flgs)
    this_flg = exp_flgs{ind_flg};    
    if flg.(this_flg)
        switch this_flg
            case 'FID_Rb'
                func = 'Alkali_FID_measurement';
                label = 'Alkali FID';
            case 'Alkali_polarization'
                func = 'estimate_alkali_polarization_ac';
                label = 'Alkali Polarization';
            case 'FID_Xenon'
                func = 'Xenon_FID';
                label = 'Xenon FID';
            case 'measure_world_rotation'
                func = 'measure_world_rotation_new';%'measure_world_rotation_new_RealTimeSave';
                label = 'World Rotation';
            case 'measure_mgnt_sens'
                func = 'measure_mgnt_sens';
                label = 'Magnetic Responsivity';
            case 'get_data_from_cDaq'
                func = 'get_data_from_cDAQ_new';
                label = 'CDAQ Readout';
            case 'get_data_from_scope'
                func = 'get_data_from_scope';
                label = 'Scope Readout';          
            case 'find_ESR_WP'
                func = 'set_to_ESR_optimization';
                label = 'ESR WP optimization';                                
        end

        try
            logger.log([label ' measurement initialized.'])
            run(func);
            logger.log([label ' measurement done.'])
        catch e
            logger.log(['An error occured in ' e.stack(1).name ' (line '  num2str(e.stack(1).line) ')']);
            logger.log(['The message was: ' e.message])
        end
    end
end



%% Save Data
logger.log('Done')

if flg.save_data
    % save_current_working_point
    save_everything
    logger.log('Saved everything')
end

%% Finish experiment notification
reset_all = 0;  % set reset to zero
if ~flg.is_part_of_loop || prm.loop_data.ind_run == prm.loop_data.n_runs_total
    playblocking(gongObj);
end

logger.close_log
