% ========================================================================
%   This script measures the drift in measured pump and probe volrages and
%   in the estimated polarizations by iteratively running the script
%   estimate_Alkali_polarization_ac.m.
% ========================================================================

%% Closing open figures and clearing command window

close all;
clc;

%% Iteratively running estimate_Alkali_polarization_ac.m

num_iters = 20;
meas = cell(num_iters, 10);

for it=1:num_iters
    disp('-----------------');
    disp(['[+] Iteration #', num2str(it)]);
    disp('-----------------');
    
    estimate_Alkali_polarization_ac;
    
    % Retrieve data from current run
    meas{it,1} = datetime;
    meas{it,2} = Vpump_mean;
    meas{it,3} = Vpump_std;
    meas{it,4} = Vprobe_mean;
    meas{it,5} = Vprobe_std;
    meas{it,6} = b;
    meas{it,7} = delta_b;
    meas{it,8} = P0;
    meas{it,9} = PRb;
    meas{it,10} = delta_PRb;
    
    set(gcf, 'Visible', 'off'); % Suppress figure display
end

% Save results
meas_filename = ['Alkali_drift_',datestr(datetime,'yyyymmddHHMM'),'.mat'];
save(meas_filename, 'meas');

%% Analyzing data

close all;
clc;

meas = load('Alkali_drift_202306250945.mat');
meas = meas.meas;

% Voltage-to-power calibration constants (see calibrate_AOM.m)
v2p_p1 = 3.54; % mW/V
v2p_p2 = 2.31; % mW

% Extract parameters from cell array
iter_vec = 1:num_iters;
time_vec = zeros(1, numel(iter_vec));
time_vec = duration(time_vec, time_vec, time_vec);
Ppump_vec = v2p_p1 * meas{1,2} + v2p_p2;
delta_PRb_vec = zeros(1, numel(iter_vec));
PRb_vec = zeros(1, numel(iter_vec));
Vprobe_mat = zeros(numel(iter_vec), numel(Ppump_vec));
for it=1:num_iters
    time_vec(it) = meas{it,1} - meas{1,1};
    PRb_vec(it) = meas{it,9};
    delta_PRb_vec(it) = meas{it,10};
    Vprobe_mat(it,:) = sort(meas{it,4});
end

[h, m, s] = hms(time_vec);
etime_mins = 60*h + m + s/60; % Convert to minutes

% Deviation from mean probe voltage for specific pump power
idx_pump = 11;
Vprobe_rel_err = (Vprobe_mat(:,idx_pump) - mean(Vprobe_mat(:,idx_pump))) ...
                  ./ mean(Vprobe_mat(:,idx_pump)) * 100; % Percent
figure; hold on; box on;
plot(etime_mins, Vprobe_rel_err, 'LineWidth', 1);
xticks(0:10:120);
xlabel('Elapsed time (Minutes)', 'FontSize', 12);
ylabel('Deviation from mean (%)', 'FontSize', 12);

% 3D plot of pump power vs. deviation from mean probe voltage as a function of elapsed time
figure; hold on; box on;
surf(Ppump_vec, etime_mins, Vprobe_mat, 'EdgeColor', 'None');
ylim([0,15]);
xticks(0:5:30);
yticks(0:5:15);
xlabel('Pump intensity (mW)', 'FontSize', 12);
ylabel('Elapsed time (Minutes)', 'FontSize', 12);
zlabel('Deviation from mean (V)', 'FontSize', 12);
view(3);

% Esimated Rb polarization as a function of elapsed time
figure; hold on; box on;
errorbar(etime_mins, PRb_vec * 100, delta_PRb_vec * 100, '.', 'LineWidth', 1);
title(['Mean Polarization: ', num2str(mean(PRb_vec)*100,4), '%,  Std: ', ...
       num2str(std(PRb_vec)*100,3), '%'], 'FontSize', 12);

xlim([0,15]);
xlabel('Elapsed time (Minutes)', 'FontSize', 12);
ylabel('Rb Polarization @ 20mW (%)', 'FontSize', 12);
