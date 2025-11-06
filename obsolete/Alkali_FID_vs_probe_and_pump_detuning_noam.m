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
%
%   TODO:
%   Add already computed errors into final plots as error bars.
%
% ========================================================================



logger.log('=======  initializing an FID rubidium vs both Lasers detuning experiment  =======')


% set pump detuning limits
det_lo_lim_pump = prm.WP.lasers.pump_detuning + prm.pump_span_volts(1);
det_up_lim_pump = prm.WP.lasers.pump_detuning + prm.pump_span_volts(2);

% set pump detuning values
meas.laser_det_pump = linspace(det_lo_lim_pump, det_up_lim_pump, prm.pump_span_points);
mm = length(meas.laser_det_pump);

% set pump detuning limits
det_lo_lim_probe = prm.WP.lasers.probe_detuning + prm.probe_span_volts(1);
det_up_lim_probe = prm.WP.lasers.probe_detuning + prm.probe_span_volts(2);

% set pump detuning values
meas.laser_det_probe = linspace(det_lo_lim_probe, det_up_lim_probe, prm.probe_span_points);
nn = length(meas.laser_det_probe);

% total number of measurements
ss = nn * mm;

% setting 2D vectors for results (mm, nn) = (pump_det, probe_det)
meas.Rb_T2_85_2d = zeros(mm, nn);
meas.Rb_T2_85_err_2d = zeros(mm, nn);
meas.Rb_a85_2d = zeros(mm, nn);
meas.Rb_a85_err_2d = zeros(mm, nn);
meas.Rb_f85_2d = zeros(mm, nn);
meas.Rb_f85_err_2d = zeros(mm, nn);
meas.Rb_phi85_2d = zeros(mm, nn);
meas.Rb_phi85_err_2d = zeros(mm, nn);

% set empty vectors for the 87 results
meas.Rb_T2_87_2d = zeros(mm, nn);
meas.Rb_T2_87_err_2d = zeros(mm, nn);
meas.Rb_a87_2d = zeros(mm, nn);
meas.Rb_a87_err_2d = zeros(mm, nn);
meas.Rb_f87_2d = zeros(mm, nn);
meas.Rb_f87_err_2d = zeros(mm, nn);
meas.Rb_phi87_2d = zeros(mm, nn);
meas.Rb_phi87_err_2d = zeros(mm, nn);

% set empty vector for another parameters results
meas.Rb_c_2d = zeros(mm, nn);
meas.Rb_c_err_2d = zeros(mm, nn);

% set probe sweep measurement
flg.FID_Rb_vs_Laser_wavelength = 1;
%chl=2;
meas.raw_data{1,1}=0;
FIDindex=1;
for i2 = 1:nn
    % set probe parameters to instrument
    AG4.DC(2, meas.laser_det_probe(i2))

    % run Alkali FID vs probe detuning
    %Alkali_FID_vs_laser_detuning
    Alkali_FID_vs_laser_detuning_noam


    % append the single iteration 85 results to matrix results 
    meas.Rb_T2_85_2d(:, i2) = meas.Rb_T2_85;
    meas.Rb_T2_85_err_2d(:, i2) = meas.Rb_T2_85_err;
    meas.Rb_a85_2d(:, i2) = meas.Rb_a85;
    meas.Rb_a85_err_2d(:, i2) = meas.Rb_a85_err;
    meas.Rb_f85_2d(:, i2) = meas.Rb_f85;
    meas.Rb_f85_err_2d(:, i2) = meas.Rb_f85_err;
    meas.Rb_phi85_2d(:, i2) = meas.Rb_phi85;
    meas.Rb_phi85_err_2d(:, i2) = meas.Rb_phi85_err;

    % append the single iteration 87 results to matrix results
    meas.Rb_T2_87_2d(:, i2) = meas.Rb_T2_87;
    meas.Rb_T2_87_err_2d(:, i2) = meas.Rb_T2_87_err;
    meas.Rb_a87_2d(:, i2) = meas.Rb_a87;
    meas.Rb_a87_err_2d(:, i2) = meas.Rb_a87_err;
    meas.Rb_f87_2d(:, i2) = meas.Rb_f87;
    meas.Rb_f87_err_2d(:, i2) = meas.Rb_f87_err;
    meas.Rb_phi87_2d(:, i2) = meas.Rb_phi87;
    meas.Rb_phi87_err_2d(:, i2) = meas.Rb_phi87_err;
    
    % append the single iteration another parameters results to matrix results
    meas.Rb_c_2d(:, i2) = meas.Rb_c;
    meas.Rb_c_err_2d(:, i2) = meas.Rb_c_err;
end

% create a meshgrid for lasers detuninigs
if 0
%[X, Y] = meshgrid([meas.laser_det_pump,meas.laser_det_pump(length(meas.laser_det_pump):-1:1)], meas.laser_det_probe);
[X, Y] = meshgrid(meas.laser_det_pump, meas.laser_det_probe);

% plotting 2D results
figure
mesh(X, Y, meas.Rb_T2_85_2d')
xlabel('Pump det (volts)')
ylabel('Probe det (volts)')
zlabel('Rb_{85} T2 (sec)')
figure
mesh(X, Y, meas.Rb_a85_2d')
xlabel('Pump det (volts)')
ylabel('Probe det (volts)')
zlabel('Rb_{85} a (volts)')
figure
mesh(X, Y, meas.Rb_f85_2d')
xlabel('Pump det (volts)')
ylabel('Probe det (volts)')
zlabel('Rb_{85} f (sec)')
figure
mesh(X, Y, meas.Rb_phi85_2d')
xlabel('Pump det (volts)')
ylabel('Probe det (volts)')
zlabel('Rb_{85} phi (rad)')
figure
mesh(X, Y, meas.Rb_T2_87_2d')
xlabel('Pump det (volts)')
ylabel('Probe det (volts)')
zlabel('Rb_{87} T2 (sec)')
figure
mesh(X, Y, meas.Rb_a87_2d')
xlabel('Pump det (volts)')
ylabel('Probe det (volts)')
zlabel('Rb_{87} a (volts)')
figure
mesh(X, Y, meas.Rb_f87_2d')
xlabel('Pump det (volts)')
ylabel('Probe det (volts)')
zlabel('Rb_{87} f (Hez)')
figure
mesh(X, Y, meas.Rb_phi87_2d')
xlabel('Pump det (volts)')
ylabel('Probe det (volts)')
zlabel('Rb_{87} phi (rad)')
end