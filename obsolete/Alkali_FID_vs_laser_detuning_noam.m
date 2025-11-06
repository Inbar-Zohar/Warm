% ========================================================================
%   This script is used for Alkali FID measurement vs laser detunings.
%   In order to choose between the pump / probe , set the experiment
%   flag flg.FID_Rb_vs_Laser_wavelength to 1 / 2 respectively.
%
%   This script uses the Alkali_FID_measurement script, and based on the 
%   assumptions written in the headline of the that script.
% ========================================================================
%
%   TODO:
%   Add error bars to resulted plots, the errors are already computed, so
%   they need only to be ploted with the data.
% ========================================================================



logger.log('=======  initializing an FID rubidium vs Laser detuning experiment  =======')

% which laser
if flg.FID_Rb_vs_Laser_wavelength == 1
    prm.which_laser = 'Pump';
    prm.laser_span_volts = prm.pump_span_volts;
    prm.laser_span_points = prm.pump_span_points;
    det_lo_lim = prm.WP.lasers.pump_detuning + prm.laser_span_volts(1);
    det_up_lim = prm.WP.lasers.pump_detuning + prm.laser_span_volts(2);
end
if flg.FID_Rb_vs_Laser_wavelength == 2
    prm.which_laser = 'Probe';
    prm.laser_span_volts = prm.probe_span_volts;
    prm.laser_span_points = prm.probe_span_points;
    det_lo_lim = prm.WP.lasers.probe_detuning + prm.laser_span_volts(1);
    det_up_lim = prm.WP.lasers.probe_detuning + prm.laser_span_volts(2);
end

% set laser detuning values
meas.laser_det = linspace(det_lo_lim, det_up_lim, prm.laser_span_points);
nnn = length(meas.laser_det);

% set empty vectors for the 85 results
meas.Rb_T2_85 = zeros(1, nnn);
meas.Rb_T2_85_err = zeros(1, nnn);
meas.Rb_a85 = zeros(1, nnn);
meas.Rb_a85_err = zeros(1, nnn);
meas.Rb_f85 = zeros(1, nnn);
meas.Rb_f85_err = zeros(1, nnn);
meas.Rb_phi85 = zeros(1, nnn);
meas.Rb_phi85_err = zeros(1, nnn);

% set empty vectors for the 87 results
meas.Rb_T2_87 = zeros(1, nnn);
meas.Rb_T2_87_err = zeros(1, nnn);
meas.Rb_a87 = zeros(1, nnn);
meas.Rb_a87_err = zeros(1, nnn);
meas.Rb_f87 = zeros(1, nnn);
meas.Rb_f87_err = zeros(1, nnn);
meas.Rb_phi87 = zeros(1, nnn);
meas.Rb_phi87_err = zeros(1, nnn);

% set empty vectors for the another parameters results
meas.Rb_c=zeros(1,nnn);
meas.Rb_c_err=zeros(1,nnn);

% loop experiments
nnn=mm;

for i1 = 1:nnn
    % set the laser detuning
    chl = flg.FID_Rb_vs_Laser_wavelength;
    AG4.DC(chl, meas.laser_det(i1))

    % puase until system stabilize
    pause(4);

    % display and log experiment
    if flg.FID_Rb_vs_Pump_And_Probe_detuning
        disp(['| meas #: (' num2str((i2 - 1) * nnn + i1) ' / ' num2str(ss) ') | Set ' prm.which_laser ' laser detuning to ' num2str(meas.laser_det(i1)) ' volts. |'])
        logger.log(['| meas #: (' num2str((i2 - 1) * nnn + i1) ' / ' num2str(ss) ') | Set ' prm.which_laser ' laser detuning to ' num2str(meas.laser_det(i1)) ' volts. |'])
    else
        disp(['| meas #: (' num2str(i1) ' / ' num2str(nnn) ') | Set ' prm.which_laser ' laser detuning to ' num2str(meas.laser_det(i1)) ' volts. |'])
        logger.log(['| meas #: (' num2str(i1) ' / ' num2str(nnn) ') | Set ' prm.which_laser ' laser detuning to ' num2str(meas.laser_det(i1)) ' volts. |'])
    end

    % run a single Alkali FID experiment
    try
        Alkali_FID_measurement_noam2;
    catch
        disp('A problem occured in fitting process.');
        logger.log('A problem occured in fitting process, skipping to next measurement.')
        continue
    end
    
    % closing all single realization plots
    if flg.close_mid_plots
       close all 
    end

    % save data
    meas.Rb_T2_85(i1) = meas.FID_T2_85;
    meas.Rb_T2_85_err(i1) = meas.FID_T2_85_err;
    meas.Rb_a85(i1) = meas.FID_a85;
    meas.Rb_a85_err(i1) = meas.FID_a85_err;
    meas.Rb_f85(i1) = meas.FID_f85;
    meas.Rb_f85_err(i1) = meas.FID_f85_err;
    meas.Rb_phi85(i1) = meas.FID_phi85;
    meas.Rb_phi85_err(i1) = meas.FID_phi85_err;
    
    meas.Rb_T2_87(i1) = meas.FID_T2_87;
    meas.Rb_T2_87_err(i1) = meas.FID_T2_87_err;
    meas.Rb_a87(i1) = meas.FID_a87;
    meas.Rb_a87_err(i1) = meas.FID_a87_err;
    meas.Rb_f87(i1) = meas.FID_f87;
    meas.Rb_f87_err(i1) = meas.FID_f87_err;
    meas.Rb_phi87(i1) = meas.FID_phi87;
    meas.Rb_phi87_err(i1) = meas.FID_phi87_err;

    meas.Rb_c(i1) = meas.FID_c;
    meas.Rb_c_err(i1) = meas.FID_c_err;
end
if 0
for i1 = nnn:-1:1
    % set the laser detuning
    chl = flg.FID_Rb_vs_Laser_wavelength;
    AG4.DC(chl, meas.laser_det(i1))

    % puase until system stabilize
    pause(4);

    % display and log experiment
    if flg.FID_Rb_vs_Pump_And_Probe_detuning
        disp(['| meas #: (' num2str((i2 - 1) * (2 * nnn+1) + i1) ' / ' num2str(2*ss) ') | Set ' prm.which_laser ' laser detuning to ' num2str(meas.laser_det(i1)) ' volts. |'])
        logger.log(['| meas #: (' num2str((i2 - 1) * (2 * nnn+1) + i1) ' / ' num2str(ss) ') | Set ' prm.which_laser ' laser detuning to ' num2str(meas.laser_det(i1)) ' volts. |'])
    else
        disp(['| meas #: (' num2str(i1) ' / ' num2str(nnn) ') | Set ' prm.which_laser ' laser detuning to ' num2str(meas.laser_det(i1)) ' volts. |'])
        logger.log(['| meas #: (' num2str(i1) ' / ' num2str(nnn) ') | Set ' prm.which_laser ' laser detuning to ' num2str(meas.laser_det(i1)) ' volts. |'])
    end

    % run a single Alkali FID experiment
    try
        Alkali_FID_measurement_noam2;
    catch
        disp('A problem occured in fitting process.');
        logger.log('A problem occured in fitting process, skipping to next measurement.')
        continue
    end
    
    % closing all single realization plots
    if flg.close_mid_plots
       close all 
    end

    % save data
    meas.Rb_T2_85(i1+nnn) = meas.FID_T2_85;
    meas.Rb_T2_85_err(i1+nnn) = meas.FID_T2_85_err;
    meas.Rb_a85(i1+nnn) = meas.FID_a85;
    meas.Rb_a85_err(i1+nnn) = meas.FID_a85_err;
    meas.Rb_f85(i1+nnn) = meas.FID_f85;
    meas.Rb_f85_err(i1+nnn) = meas.FID_f85_err;
    meas.Rb_phi85(i1+nnn) = meas.FID_phi85;
    meas.Rb_phi85_err(i1+nnn) = meas.FID_phi85_err;
    
    meas.Rb_T2_87(i1+nnn) = meas.FID_T2_87;
    meas.Rb_T2_87_err(i1+nnn) = meas.FID_T2_87_err;
    meas.Rb_a87(i1+nnn) = meas.FID_a87;
    meas.Rb_a87_err(i1+nnn) = meas.FID_a87_err;
    meas.Rb_f87(i1+nnn) = meas.FID_f87;
    meas.Rb_f87_err(i1+nnn) = meas.FID_f87_err;
    meas.Rb_phi87(i1+nnn) = meas.FID_phi87;
    meas.Rb_phi87_err(i1+nnn) = meas.FID_phi87_err;

    
end
end
% return to initilaized laser detuning value
if chl == 1
    AG4.DC(chl, prm.WP.lasers.pump_detuning);
end
if chl == 2
    AG4.DC(chl, prm.WP.lasers.probe_detuning);
end
flg.plot_data=1;

if flg.plot_data

    figure;
    subplot(211)
    plot(meas.laser_det_pump, meas.Rb_T2_85(1:nnn), 'r.', 'markerSize', 8);
     hold on;
     plot(meas.laser_det_pump, meas.Rb_T2_87(1:nnn), 'b.', 'markerSize', 8);
     hold off;
    legend('T_2', 'Location', 'Best');
    xlabel('Laser detuning (volts)')
    ylabel('T_2 (sec)')
    title(['Rubidium 85 T_2 vs pump laser detuning with probe laser detuning=' num2str(meas.laser_det_probe(i2))])
    subplot(212)
    plot(meas.laser_det_pump, meas.Rb_a85(1:nnn), 'r.', 'markerSize', 10);
     hold on;
     plot(meas.laser_det_pump, meas.Rb_a87(1:nnn), 'b.', 'markerSize', 10);
     hold off;
    title(['Rubidium 85 amplitude vs pump laser detuning with probe laser detuning=' num2str(meas.laser_det_probe(i2))])
    xlabel('Laser detuning (volts)')
    ylabel('Amplitude (volts)')
    legend('Amplitude', 'Location', 'Best');
end

logger.log('=======  FID rubidium vs Laser wavelength experiment is done.  =======')
