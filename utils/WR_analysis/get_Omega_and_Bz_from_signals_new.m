function [Omega_deg_hr, dBz, Omega_129_deg_hr, Omega_131_deg_hr,dESR] = ...
    get_Omega_and_Bz_from_signals_new(meas,prm)

%% Contsants
ga = (-7.44e7  )/(2*pi)/(1e4);              % [Hz/G] Xe129 gyromagnetic ratio
gb = (+2.2056e7)/(2*pi)/(1e4);              % [Hz/G] Xe131 gyromagnetic ratio
x = ga / gb;
R = abs(x);

%% Data preparation and legacy shit
if isfield(prm,'locking_Bz_instead_of_drive')
    if sum(prm.locking_Bz_instead_of_drive) == 2
        warning('Error: You can''t have both PIDs locking Bz. Didn''t calculte Omega')
        Omega_deg_hr = nan(size(meas.time_s));
        dBz = nan(size(meas.time_s));
        Omega_129_deg_hr = nan(size(meas.time_s));
        Omega_131_deg_hr = nan(size(meas.time_s));
    else 
        locking_Bz_instead_of_drive = prm.locking_Bz_instead_of_drive;
    end
elseif isfield(prm,'american_lock') && prm.american_lock
    locking_Bz_instead_of_drive = [0 0];
else
    locking_Bz_instead_of_drive = [1 0];
end
no_131_signal = 0;
if isfield(meas,'PID_131_out')
    PID_131_out = meas.PID_131_out;
elseif isfield(meas,'Det131_PIDb_V')
    PID_131_out = meas.Det131_PIDb_V;
else
    no_131_signal = 1;
    PID_131_out = nan(size(meas.time_s));
end

no_129_signal = 0;
if isfield(meas,'PID_129_out')
    PID_129_out = meas.PID_129_out;
elseif isfield(meas,'Bz_V')
    PID_129_out = meas.Bz_V;
else
    no_129_signal = 1;    
    PID_129_out = nan(size(meas.time_s));
end

%% Actual calculation
if no_131_signal || no_129_signal
    %% Missing channel
    if no_131_signal && ~no_129_signal
        warning('Error: Missing the 131 PID output channel. Didn''t calculte Omega')
    elseif ~no_131_signal && no_129_signal
        warning('Error: Missing the 129 PID output channel. Didn''t calculte Omega')
    elseif no_131_signal && no_129_signal
        warning('Error: Missing both 129,131 PID output channels. Didn''t calculte Omega')
    end
else
    if ~any(locking_Bz_instead_of_drive)
        %% Locking both NMR signals
        Omega_129_Hz = prm.WP.Xe129.mod_dev / 5 * PID_129_out;
        Omega_131_Hz = prm.WP.Xe131.mod_dev / 5 * PID_131_out;

        
        Omega_129_deg_hr = Omega_129_Hz * 3600*360;
        Omega_131_deg_hr = Omega_131_Hz * 3600*360;

        Omega_deg_hr = Omega_129_deg_hr/(1+R)-Omega_131_deg_hr*R/(1+R);
        
        dESR = (Omega_129_Hz+Omega_131_Hz)/(abs(ga)+gb);  % [G]

        Bz_calib =  0.0316 * 1e-3 * prm.calib.coils.dBz_dIz;%/(abs(ga)+gb); %  30dB attenuator * 1mA/V * 2.98G/A
        dBz = Bz_calib * meas.PID_simR_out; 
    else
        %% Locking one NMR and Bz
        % Bz_calib = -100e-6*40/150*3; % CHECK CALIBRATION 100uA/V*40uA/150uA*3G/A!
        PID_sim_out = meas.PID_sim_out;
        Bz_calib =  0.0316 * 1e-3 * prm.calib.coils.dBz_dIz;%/(abs(ga)+gb); %  30dB attenuator * 1mA/V * 2.98G/A 
        Omega_calib =  1 / 5 ;%/ (1 - 1 / x);
        if locking_Bz_instead_of_drive(1)
            Bz_sig = PID_sim_out;
            Omega_sig131 = prm.WP.Xe131.mod_dev*PID_131_out;
            Omega_sig129 = prm.WP.Xe129.mod_dev*PID_129_out;
        else
            Bz_sig = PID_131_out;
            Omega_sig = PID_129_out;
        end
        dBz = Bz_calib * Bz_sig;% [G] 
        Omega_deg_hr131 =  (abs(gb)/(2*pi) * dBz+Omega_calib *Omega_sig131) * 360 * 3600;%                               
        Omega_deg_hr129 =  (abs(ga)/(2*pi) * dBz+Omega_calib *Omega_sig129) * 360 * 3600;%                               
         
        
        Omega_129_deg_hr = Omega_deg_hr129;% abs(ga)/(2*pi) * dBz + Omega_deg_hr129;
        Omega_131_deg_hr = Omega_deg_hr131;%abs(gb)/(2*pi) * dBz - Omega_deg_hr131;
        Omega_deg_hr = Omega_129_deg_hr/(1+R)-Omega_131_deg_hr*R/(1+R);
        % Omega_deg_hr = Omega_deg_hr131;
        warning('Check Omega,Bz->Omega129,131 signs!');
    end
end