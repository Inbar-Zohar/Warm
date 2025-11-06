% Automation wrapper for: main_.m
%
% Note: limited functionality available at the moment!
%
% Before running this file (which runs main_.m multiple times), all the
% configurations in main_.m should be set for the experiment we want to run.
%
%
% Multiple measurements configurations
%
% Syntax:
% mult - cell array. Each cell represents a seperate experiment
%
% mult{I}.dev_list - cell array of the devices we want to configure in an
% experiment
%
% mult{I}.dev_list{J} - struct with the following fields:
%   dev:        A reference to the device to configure. (For example
%               mult{1}.dev_list{3}.dev = AG1).
%               At the moment only AG1-6 are supported (Kesight 33500B series)
%   dev_name:   The device name, for use in the automatic experiment
%               description. (For example: mult{1}.dev_list{3}.dev_name =
%               'AG1')
%   dev_ch:     The channel to configure
%   dev_waveform:   The waveform. At the moment only 'DC','SIN','RAMP' are
%                   supported
%   dev_par_vec:    The parameters vector. Each waveform should have a
%                   parameters vector according to its syntax:
%                   'DC': bias (V)
%                   'SIN': [freq. (Hz), amp. (Vpp), offset (V)]
%                   'RAMP': [freq. (Hz), amp. (Vpp), offset (V), phase (deg)]
%   
%
% mult_no_run_indices: A ROW vector of experiment indices for which we don't
%                      want to run 'main_.m'. This is useful in cases where
%                      we want to change device parameters but not run the
%                      experiment. For example: useful to set the
%                      configurations back at the end of all experiments.
%



%% Initializations
% if exist('mult','var')
    %Make sure to delete 'mult' so we don't accidently use configurations
    %of a previous run.
    clearvars mult mult_template mult_no_run_indices;
% end

%% Configurations

%What wait time to use, to let the system stabalize between measurements:
WAIT_TIME = 60; %sec

% mult_all_ramp_amps = [logspace(-1.15,0.45,8) 5]; % Vpp
% mult_all_ramp_amps = (10:10:80)*1e-3;
mult_all_ramp_amps = [linspace(8,80,10)*1e-3 ... % experiments set #1 - Gy
                      linspace(1,8,10)*1e-3 ... % experiments set #2 - Pump
                      linspace(1,8,10)*1e-3 ... % experiments set #3 - Probe
                     ];
% exp_order = [1 2 1 2 1 3 1 3];
exp_order = ones([1 60]);
exp_order(2:2:20) = 2;
exp_order(22:2:40) = 3;
exp_order(42:2:60) = 4;

%Exp 1 (template for all experiments without a ramp)
Ie = 1;
Id = 1;
mult_template{Ie}.dev_list{Id}.dev = AG4; % Reference to a device
mult_template{Ie}.dev_list{Id}.dev_name = 'AG4'; % Name of the device
mult_template{Ie}.dev_list{Id}.dev_ch = 1; % Device channel number
mult_template{Ie}.dev_list{Id}.dev_waveform = 'DC'; % The waveform
mult_template{Ie}.dev_list{Id}.dev_par_vec = 175e-3; % offset (V)
Id = 2;
mult_template{Ie}.dev_list{Id}.dev = AG4; % Reference to a device
mult_template{Ie}.dev_list{Id}.dev_name = 'AG4'; % Name of the device
mult_template{Ie}.dev_list{Id}.dev_ch = 2; % Device channel number
mult_template{Ie}.dev_list{Id}.dev_waveform = 'DC'; % The waveform
mult_template{Ie}.dev_list{Id}.dev_par_vec = -2e-3; % offset (V)
Id = 3;
mult_template{Ie}.dev_list{Id}.dev = AG3; % Reference to a device
mult_template{Ie}.dev_list{Id}.dev_name = 'AG3'; % Name of the device
mult_template{Ie}.dev_list{Id}.dev_ch = 2; % Device channel number
mult_template{Ie}.dev_list{Id}.dev_waveform = 'DC'; % The waveform
mult_template{Ie}.dev_list{Id}.dev_par_vec = 0; % offset (V)

%Exp 2 (template for all experiments with a ramp)
Ie = 2;
Id = 1;
mult_template{Ie}.dev_list{Id}.dev = AG3; % Reference to a device
mult_template{Ie}.dev_list{Id}.dev_name = 'AG3'; % Name of the device
mult_template{Ie}.dev_list{Id}.dev_ch = 2; % Device channel number
mult_template{Ie}.dev_list{Id}.dev_waveform = 'RAMP'; % The waveform
mult_template{Ie}.dev_list{Id}.dev_par_vec = [0.272e-3 1e-3 0 180]; % 

%Exp 3 (template for all experiments with a ramp)
Ie = 3;
Id = 1;
mult_template{Ie}.dev_list{Id}.dev = AG4; % Reference to a device
mult_template{Ie}.dev_list{Id}.dev_name = 'AG4'; % Name of the device
mult_template{Ie}.dev_list{Id}.dev_ch = 1; % Device channel number
mult_template{Ie}.dev_list{Id}.dev_waveform = 'RAMP'; % The waveform
mult_template{Ie}.dev_list{Id}.dev_par_vec = [0.272e-3 1e-3 175e-3 180]; % 

%Exp 4 (template for all experiments with a ramp)
Ie = 4;
Id = 1;
mult_template{Ie}.dev_list{Id}.dev = AG4; % Reference to a device
mult_template{Ie}.dev_list{Id}.dev_name = 'AG4'; % Name of the device
mult_template{Ie}.dev_list{Id}.dev_ch = 2; % Device channel number
mult_template{Ie}.dev_list{Id}.dev_waveform = 'RAMP'; % The waveform
mult_template{Ie}.dev_list{Id}.dev_par_vec = [0.272e-3 1e-3 -2e-3 180]; % 


for Ic = 1:numel(mult_all_ramp_amps)
    mult{2*Ic-1} = mult_template{exp_order(2*Ic-1)}; %DC experiment
    mult{2*Ic} = mult_template{exp_order(2*Ic)}; % Ramp experiment
    mult{2*Ic}.dev_list{Id}.dev_par_vec(2) = mult_all_ramp_amps(Ic); %Update value
end

%At the end set the device back to DC:
mult{2*Ic+1} = mult_template{1}; %DC experiment
mult_no_run_indices = [2*Ic+1];
% mult_no_run_indices = 1:numel(mult);

%% Script
try

N_skip = sum(sum( repmat( (1:numel(mult)), [numel(mult_no_run_indices) 1] ) == mult_no_run_indices.' ));
N_exp = numel(mult) - N_skip;
I_exp = 1;

for Imult = 1:numel(mult)
 %Iterate over all experiments
%  disp(num2str(Imult));
 
 %Create an automatic experiment description
 mult_cur_exp_description = ['; ' num2str(I_exp) '/' num2str(N_exp) ': '];
 
 for Jmult = 1:numel( mult{Imult}.dev_list )
     %Iterate over all devices in this experiment

     %Update the automatic experiment description
     mult_cur_exp_description = [mult_cur_exp_description ...
         mult{Imult}.dev_list{Jmult}.dev_name ...
         '-ch.' num2str( mult{Imult}.dev_list{Jmult}.dev_ch ) ...
         ' - ' mult{Imult}.dev_list{Jmult}.dev_waveform]; %#ok<*AGROW>
     
     %Set the device parameters (waveform)
     switch lower( mult{Imult}.dev_list{Jmult}.dev_waveform )
         case 'dc'
             mult{Imult}.dev_list{Jmult}.dev.DC(...
                 mult{Imult}.dev_list{Jmult}.dev_ch, ... %channel
                 mult{Imult}.dev_list{Jmult}.dev_par_vec(1) ... % offset (V)
             );
            %Automatically change the experiment description
            mult_cur_exp_description = [mult_cur_exp_description ...
                   ', ' num2str( mult{Imult}.dev_list{Jmult}.dev_par_vec(1) ) ' V' ' | '];

         case 'sin'
            mult{Imult}.dev_list{Jmult}.dev.Sin(...
                mult{Imult}.dev_list{Jmult}.dev_ch, ... %channel
                mult{Imult}.dev_list{Jmult}.dev_par_vec(1), ... % freq. (Hz)
                mult{Imult}.dev_list{Jmult}.dev_par_vec(2), ... % amp. (Vpp)
                0, ... % phase (deg) - hardcoded.
                mult{Imult}.dev_list{Jmult}.dev_par_vec(3) ... % offset (V)
            );
            %Automatically change the experiment description
            mult_cur_exp_description = [mult_cur_exp_description ...
                   ', ' num2str( mult{Imult}.dev_list{Jmult}.dev_par_vec(1) ) ' Hz, ' ...
                   num2str( mult{Imult}.dev_list{Jmult}.dev_par_vec(2) ) ' Vpp, ' ...
                   num2str( mult{Imult}.dev_list{Jmult}.dev_par_vec(3) ) ' V' ' | '];

         case 'ramp'
            mult{Imult}.dev_list{Jmult}.dev.Ramp(...
                mult{Imult}.dev_list{Jmult}.dev_ch, ... %channel
                mult{Imult}.dev_list{Jmult}.dev_par_vec(1), ... % freq. (Hz)
                mult{Imult}.dev_list{Jmult}.dev_par_vec(2), ... % amp. (Vpp)
                mult{Imult}.dev_list{Jmult}.dev_par_vec(3), ... % offset (V)
                mult{Imult}.dev_list{Jmult}.dev_par_vec(4) ... % phase (deg).
            );
            %Automatically change the experiment description
            mult_cur_exp_description = [mult_cur_exp_description ...
                   ', ' num2str( mult{Imult}.dev_list{Jmult}.dev_par_vec(1) ) ' Hz, ' ...
                   num2str( mult{Imult}.dev_list{Jmult}.dev_par_vec(2) ) ' Vpp, ' ...
                   num2str( mult{Imult}.dev_list{Jmult}.dev_par_vec(3) ) ' V, ' ...
                   num2str( mult{Imult}.dev_list{Jmult}.dev_par_vec(4) ) ' Deg' ' | '];

         otherwise
            %TODO: Add support for other waveforms as required!
            disp(['The ' mult{Imult}.dev_list{Jmult}.dev_waveform ' waveform is not supported in multiple_measurements mode.'])
            logger.log(['The ' mult{Imult}.dev_list{Jmult}.dev_waveform ' waveform is not supported in multiple_measurements mode.'])                            
     end
     
 
     
 end
 
 %Allow the system a few seconds to "settle" before starting the measuerment
 pause(WAIT_TIME);

 %Run the experiment
 mult_flg_is_running = 1;
 close all;
 
 if ~sum(mult_no_run_indices == Imult)
    I_exp = I_exp + 1;
%     disp(['run' mult_cur_exp_description]);
    run('main_.m');
 end
 
 mult_flg_is_running = 0;
 
 
end
    
catch e
    disp(['a problem occured in multiple_measurements mode. The message was: ' e.message])
    logger.log(['a problem occured in multiple_measurements mode. The message was: ' e.message]) 
end

