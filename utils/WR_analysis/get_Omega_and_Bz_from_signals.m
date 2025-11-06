function [Omega_deg_hr, dBz, Omega_129_deg_hr, Omega_131_deg_hr] = ...
    get_Omega_and_Bz_from_signals(meas,prm)

ga = (-7.44e7  )/(2*pi)/(1e4);              % [Hz/G] Xe129 gyromagnetic ratio
gb = (+2.2056e7)/(2*pi)/(1e4);              % [Hz/G] Xe131 gyromagnetic ratio
x = ga / gb;

% data

if ~isfield(prm,'american_lock') || ~prm.american_lock    
    W_Hz = prm.df131 / 5 / (1 - 1 / x) * meas.Det131_PIDb_V;      % [Hz]
    Omega_deg_hr = W_Hz * 360 * 3600;                               % [deg/hr]
    dBz = -100e-6*40/150*3 * meas.Bz_V;                               % [G] 100uA/V*40uA/150uA*3G/A!      
    Omega_129_deg_hr = nan(size(Omega_deg_hr));
    Omega_131_deg_hr = nan(size(Omega_deg_hr));
else  
    % american lock from improvement of stability
    %     warning('american lock now')
    R = abs(x);
    W_Hz_131 = prm.df131 / 5 * meas.Det131_PIDb_V;      % [Hz]
    W_Hz_129 = prm.df129 / 5 * meas.Bz_V;      % [Hz]
    W_Hz = W_Hz_129/(1+R)-W_Hz_131*R/(1+R);
    
    Omega_129_deg_hr = W_Hz_129 * 3600*360;
    Omega_131_deg_hr = W_Hz_131 * 3600*360;
    
    Omega_deg_hr = W_Hz * 360 * 3600;
    dBz = (W_Hz_129+W_Hz_131)/(abs(ga)+gb);  % [G] 100uA/V*40uA/150uA*3G/A!    
end
