function meas = extract_raw_data_from_MFLI(file)

tmp = load(file);
dev32531 = tmp.dev32531;
meas.Xe_131_X = dev32531.demods(3).sample.x;
meas.Xe_129_X = dev32531.demods(4).sample.x;

meas.Xe_131_Y = dev32531.demods(3).sample.y;
meas.Xe_129_Y = dev32531.demods(4).sample.y;

meas.Xe_131_Theta = atan2(meas.Xe_131_X,meas.Xe_131_Y) * 180 / pi;
meas.Xe_129_Theta = atan2(meas.Xe_129_X,meas.Xe_129_Y) * 180 / pi;

meas.Xe_131_R = sqrt(meas.Xe_131_X.^2+meas.Xe_131_Y.^2);
meas.Xe_129_R = sqrt(meas.Xe_129_X.^2+meas.Xe_129_Y.^2);

Omega_129_Hz = dev32531.demods(3).sample.frequency - mean(dev32531.demods(3).sample.frequency) ;
Omega_131_Hz = dev32531.demods(4).sample.frequency - mean(dev32531.demods(4).sample.frequency) ;

meas.Omega_129_deg_hr = Omega_129_Hz * 360 * 3600;
meas.Omega_131_deg_hr = Omega_131_Hz * 360 * 3600;



ga = (-7.44e7  )/(2*pi)/(1e4);              % [Hz/G] Xe129 gyromagnetic ratio
gb = (+2.2056e7)/(2*pi)/(1e4);              % [Hz/G] Xe131 gyromagnetic ratio
x = ga / gb;
R = abs(x);

meas.Omega_deg_hr = meas.Omega_129_deg_hr/(1+R)-meas.Omega_131_deg_hr*R/(1+R);
meas.Bz = (Omega_129_Hz+Omega_131_Hz)/(abs(ga)+gb);  % [G]

meas.time_s = double((dev32531.demods(3).sample.timestamp - dev32531.demods(3).sample.timestamp(1)))/60e6;
meas.fs = 1./diff(meas.time_s(1:2));
[~,timestamp_str] = fileparts(fileparts(fileparts(file)));
[meas.DateTime_var, meas.DateNumber] = extract_datenum_from_string(timestamp_str);


function [dt,dn] = extract_datenum_from_string(str)
    % Extract the relevant date-time portion from the string
    % Example: 'session_20250817_100115_33'
    
    % Split the string by underscores
    parts = split(str, '_');
    
    % Extract date and time parts
    datePart = parts{2}; % '20250817'
    timePart = parts{3}; % '100115'
    
    % Combine into datetime string
    dtStr = [datePart timePart]; % '20250817100115'
    
    % Convert to datetime
    dt = datetime(dtStr, 'InputFormat', 'yyyyMMddHHmmss');
    
    % Convert to datenum
    dn = datenum(dt);
end
end