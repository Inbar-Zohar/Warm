% ========================================================================
%   This script is used for Photon Shot Noise measurement.
%
%
% ========================================================================
%
%   TODO:
%   measure the optical power directly from the TDS2014C scope from the
%   monitor channels.
% ========================================================================

logger.log('LIA noise measurement initialized.')
disp('[+] Starting LIA noise measurement');

% Configure scope to take data
chnl = 2;
Tscale = 0.2; % sec per dev
Vscale = 0.02; % volt per dev

% set scope
curr_inst = scope1;
curr_inst.Stop;

% Coupling
curr_inst.setChCoupling(chnl, 'DC');
curr_inst.setChImpedance(chnl, 'M');

% Acquire
curr_inst.HighRes;
% curr_inst.AverageAQN(2); % for 12bit data (not working as is)

% horizontal
scope1.setTref('CENTer');
curr_inst.setTscale(Tscale);
curr_inst.setTdelay(5.00 * Tscale);

% vertical offset
curr_inst.setVoffset(chnl, 0);

% vertical scale
curr_inst.setVscale(chnl, Vscale);
curr_inst.Run;

% read data from scope
curr_inst.Single;
curr_inst.readyToRead;
[t, ~, ~] = curr_inst.Read(chnl);
meas.LIA_noise_raw_time = t;
length_t = size(t, 2);
dt = mean(diff(meas.LIA_noise_raw_time));
fs = 1 / dt;

% Release scope from Single mode
curr_inst.Run; 
clear t 

% set vectors
meas.LIA_noise_raw_volts = zeros(prm.LIA_measurements, length_t);

% measure
for i = 1:prm.psn_measurements
    fprintf('Measurement num: %d:%d\n', i, prm.psn_measurements)
    % input the Probe laser optical power
    meas.shot_noise_optical_power(i) = input('What is the optical power [mW]?') * 1e-3;
    
    for j = 1:prm.psn_reps
        curr_inst.Single;
        curr_inst.readyToRead
        [~, meas.shot_noise_raw_volts(i, j, :)] = curr_inst.Read(chnl);
        curr_inst.Run;
    end
end



%% PSD calculation

% useful parameters

hc = 1.98644586e-25;                %   [J x m]
probe_wavelength = 795e-9;          %   [m]
hf = hc / probe_wavelength;         %   [J]
PDB_responsivity = 0.537;           %   [A/W]
PDB_gain = 51e3;                    %   [V/A]

% set vectors for PSD calculations
dV = squeeze(meas.shot_noise_raw_volts(1, 1, :));
[PSD, f] = pwelch(dV - mean(dV), [], [], [], fs);
PSD_size = length(PSD);
meas.PSD_data = zeros(prm.psn_measurements, prm.psn_reps, PSD_size);
meas.phi_PSD_white = zeros(prm.psn_measurements, prm.psn_reps, 1);

figure();
title('Laser power vs PSD', 'fontweight', 'bold', 'FontSize', 40, 'FontName', 'Cambria Math');
for i = 1:prm.psn_measurements
    
    % set subplot parameters
    subplot(sqrt(prm.psn_measurements), sqrt(prm.psn_measurements), i);
    for j = 1:prm.psn_reps
        
        % compute the PSD for all reps with a single optical power
        dV = squeeze(meas.shot_noise_raw_volts(i, j, :)); % RF channel 
        [meas.PSD_data(i, j, :), f] = pwelch(dV - mean(dV), [], [], [], fs);
        PSD_data_single = squeeze(meas.PSD_data(i, j, :));
        
        % plot
        loglog(f, sqrt(PSD_data_single), 'displayname', num2str(j));
        hold on
        
    end
    title(sprintf('Laser power: %.4f [W]', meas.shot_noise_optical_power(i)));
    xlabel('f(Hz)');
    ylabel('$\sqrt{PSD} (\frac{V}{\sqrt{Hz}})$', 'interpreter', 'latex');
    grid on
    legend show
    if i == 1
       % limits of white noise in spectrum
        [x_in, y_in] = ginput(2);
        white_low_f_cutoff = min(x_in);
        white_upper_f_cutoff = max(x_in); 
    end
end
close

% compute the mean of PSD over white noise spectrum
for i = 1:prm.psn_measurements
    for j = 1:prm.psn_reps        
        PSD_data_single = squeeze(meas.PSD_data(i, j, :));
        meas.phi_PSD_white(i, j) = mean(sqrt(PSD_data_single(f > ...
            white_low_f_cutoff & f < white_upper_f_cutoff)));
    end
end

% plot mosaic of the meausrements with estimated white noise bounds
figure();
title('Laser power vs PSD', 'fontweight', 'bold', 'FontSize', 40, 'FontName', 'Cambria Math');
for i = 1:prm.psn_measurements
    
    % set subplot parameters
    subplot(sqrt(prm.psn_measurements), sqrt(prm.psn_measurements), i);
    for j = 1:prm.psn_reps
        
        PSD_data_single = squeeze(meas.PSD_data(i, j, :));
        
        % plot
        loglog(f, sqrt(PSD_data_single));
        hold on
        
    end
    loglog([white_low_f_cutoff, white_upper_f_cutoff], [meas.phi_PSD_white(i, j), meas.phi_PSD_white(i, j)], 'k--')
    title(sprintf('Laser power: %.4f [W]', meas.shot_noise_optical_power(i)));
    xlabel('f(Hz)');
    ylabel('$\sqrt{PSD} (\frac{V}{\sqrt{Hz}})$', 'interpreter', 'latex');
    grid on
end


% fit the data with a Photon-Shot-Noise model
shot_noise_model = fittype('sqrt(x) * a + b', 'coefficients', ...
    {'a', 'b'});
try
    shot_noise_fit = fit(meas.shot_noise_optical_power, mean(meas.phi_PSD_white, 2), ...
        shot_noise_model, 'StartPoint', [PDB_gain * PDB_responsivity * ...
        sqrt(2 * hf), 0]);
    fit_success = 1;
catch
    disp('<8><8><8> fit did not succeed <8><8><8>')
    fit_success = 0;
end

% compute the theoretical Photon-Shot-Noise curve given the measured
% optical power
dV_theo = PDB_gain * PDB_responsivity * sqrt(2 * hf * ...
    (meas.shot_noise_optical_power));   % V * sqrt(s) = V / sqrt(Hz)

% plot the Photon-Shot-Noise curves: data vs theory
err = std(meas.phi_PSD_white, 1, 2);
figure(); hold on;
errorbar(meas.shot_noise_optical_power, mean(meas.phi_PSD_white, 2), err, 'o')
plot(meas.shot_noise_optical_power, dV_theo, 'o')

if fit_success
    plot(shot_noise_fit)
    legend('Meas', 'Theo', 'fit')
    % compute the quantum efficiency of the PDB
    BPD_quantum_eff = shot_noise_fit.a / (PDB_gain * PDB_responsivity * sqrt(2 * hf));
else
    legend('Meas', 'Theo')
end

xlabel('Probe optical power after cell [mW]', 'interpreter', 'latex')
ylabel('dV [$V/\sqrt(Hz)$]', 'interpreter', 'latex')
grid()




