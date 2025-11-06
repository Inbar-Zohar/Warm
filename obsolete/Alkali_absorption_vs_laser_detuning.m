% ========================================================================
%   This script is used for Alkali absorption measurement vs laser detunings.
%   In order to choose between the pump / probe , set the experiment
%   flag flg.Rb_absorption_vs_Laser_wavelength to 1 / 2 respectively.
%
% ========================================================================



logger.log('=======  initializing an rubidium absorption vs Laser detuning experiment  =======')

% which laser
if flg.Rb_absorption_vs_Laser_wavelength == 1
    prm.which_laser = 'Pump';
    prm.laser_span_volts = prm.pump_span_volts;
    prm.laser_span_points = prm.pump_span_points;
    det_lo_lim = prm.WP.lasers.pump_detuning + prm.laser_span_volts(1);
    det_up_lim = prm.WP.lasers.pump_detuning + prm.laser_span_volts(2);
end
if flg.Rb_absorption_vs_Laser_wavelength == 2
    prm.which_laser = 'Probe';
    prm.laser_span_volts = prm.probe_span_volts;
    prm.laser_span_points = prm.probe_span_points;
    det_lo_lim = prm.WP.lasers.probe_detuning + prm.laser_span_volts(1);
    det_up_lim = prm.WP.lasers.probe_detuning + prm.laser_span_volts(2);
end

% set laser detuning values
meas.laser_det = linspace(det_lo_lim, det_up_lim, prm.laser_span_points);
nn = length(meas.laser_det);

% set empty vectors for the results
meas.PD_diff = zeros(1, nn);
meas.PD_minus = zeros(1, nn);
meas.PD_plus = zeros(1, nn);


% loop experiments
for J = 1:nn
    % set the laser detuning
    chl = flg.Rb_absorption_vs_Laser_wavelength;
    AG4.DC(chl, meas.laser_det(J))

    % puase until system stabilize
    pause(4);

    % display and log experiment
    disp(['| meas #: (' num2str(J) ' / ' num2str(nn) ') | Set ' prm.which_laser ' laser detuning to ' num2str(meas.laser_det(J)) ' volts. |'])
    logger.log(['| meas #: (' num2str(J) ' / ' num2str(nn) ') | Set ' prm.which_laser ' laser detuning to ' num2str(meas.laser_det(J)) ' volts. |'])


    % run a single Alkali absorption experiment
    try
    % --------------------------------------------------------------------    
        logger.log('Alkali absorption measurement initialized.')
        disp('[+] Starting an Alkali absorption measurement');
        
        % Disable NMR 129 and 131 and ESR
        AG1.OutputOFF(1);   AG2.OutputOFF(1);   AG3.OutputOFF(1); AG6.OutputOFF(1);
        % We don't Disable channels 2 of AG1 and 2! the user should set it
        % to desired values.
%         AG1.OutputOFF(2); AG2.OutputOFF(2);

        % Set  AG6 ch2 to trigger all with a square pulse
        AG6.Square(2, 133, 1.5, 1.5, 0); AG6.Termination(2,50);  AG6.OutputON(2);
        
        % Set scope 1 ch1 with trigger on ch1
        % Configure scope to take data
        chn_diff = 1;
        chn_minus = 3;
        chn_plus = 4;
        Tscale = 20e-6; % sec per dev
        Vscale_initial = [0.5 1 1]; % volt per dev, for each channel
        Voffsets = [0 1.75 1.75]; % volt per dev, for each channel
        Probe_scales = [1 0.2 0.2];

        % set scope
        curr_inst = scope1;
        curr_inst.Stop;

        % Coupling
        curr_inst.setChCoupling(chn_diff, 'DC');
        curr_inst.setChImpedance(chn_diff, 'M');

        curr_inst.setChCoupling(chn_minus, 'DC');
        curr_inst.setChImpedance(chn_minus, 'M');
        
        curr_inst.setChCoupling(chn_plus, 'DC');
        curr_inst.setChImpedance(chn_plus, 'M');

        
        % Acquire
        % curr_inst.HighRes;
        % curr_inst.AverageAQN(2); % for 12bit data (not working as is)

        % horizontal
        curr_inst.setTscale(Tscale);
        % curr_inst.setTdelay(5.00 * Tscale);

        %Since the PD voltage is high, the dynamic range is limited, so for
        %the measurement we set the Probe ratio to 0.2:1 for chn_minus, and
        %chn_plus, and we will multiply the voltage measurements by 5.
        curr_inst.Probe(chn_minus, Probe_scales(2));
        curr_inst.Probe(chn_plus, Probe_scales(3));
                
        %Set initial Vscales
        curr_inst.setVscale(chn_diff, Vscale_initial(1));
        curr_inst.setVscale(chn_minus, Vscale_initial(2));
        curr_inst.setVscale(chn_plus, Vscale_initial(3));

        % vertical offsets
        curr_inst.setVoffset(chn_diff, Voffsets(1));
        curr_inst.setVoffset(chn_minus, Voffsets(2));
        curr_inst.setVoffset(chn_plus, Voffsets(3));
        

        % setting the scope vscale dynamically for each channel
        channels = [chn_diff chn_minus chn_plus];
        %To make sure absorption is minimal before setting the scale, we
        %Iterate over all channels, read data, and rescale only the plus
        %and minus channels if needed.
      for I=1:3
        
        cur_ch = channels(I);
        try_counter = 0;
        in_once_flag = 0;
        out_once_flag = 0;
        Vscale = Vscale_initial(I);

        while (try_counter < 30 && ~in_once_flag && Vscale >= 0.002)
            % setting signal to scope bounds by enlarging the scale until signal is
            % out of scope and than setting it to smaller scale until signal is in
            % again.
            
            % vertical scale
            curr_inst.setVscale(cur_ch, Vscale);
            curr_inst.Run;

            %For chn_diff we don't want to change the scale
            if (cur_ch == chn_diff)
                in_once_flag = 1;
                out_once_flag = 1;    
            else
                %For the other channels, make sure each time that the
                %offsets are set correctly
%                 curr_inst.setVoffset(cur_ch, Voffsets(I));
            end
            
            % Trigger
            curr_inst.TrigMode('NORM');
            curr_inst.TrigSource('EXT');
            curr_inst.TrigSlope('POS');
            curr_inst.TrigThresh(1);
            curr_inst.Single;

            % read data from scope
            curr_inst.readyToRead(10);
            [t, v, ~] = curr_inst.Read(cur_ch);

            in_flag = curr_inst.sigInbounds(v, cur_ch);
            if ~in_flag
                % signal not contained in scope, enlarge scale
                Vscale = Vscale * 1.75;
                try_counter = try_counter + 1;
                out_once_flag = 1;

            elseif in_flag && ~out_once_flag
                % signal is contained in scope from the start, set smaller scale
                Vscale = Vscale * 0.5;
                try_counter = try_counter + 1;

            elseif in_flag && out_once_flag
                in_once_flag = 1;
            end


        end

        % Release scope from Single mode
        curr_inst.Run;
        
        %Save the data
        if I==1
            meas.PD_diff(J) = mean(v);
        elseif I==2
            meas.PD_minus(J) = mean(v) / Probe_scales(2);
        else
            meas.PD_plus(J) = mean(v) / Probe_scales(3);
        end
        
      end
       
        
    % --------------------------------------------------------------------    
    catch
        disp('A problem occured in fitting process.');
        logger.log('A problem occured in fitting process, skipping to next measurement.')
        continue
    end
    
    % closing all single realization plots
    if flg.close_mid_plots
       close all 
    end
    
end

%The final data
meas.PD_abs = ( sqrt(meas.PD_minus) + sqrt(meas.PD_plus) ).^2;
% meas.PD_abs = ( sqrt(meas.PD_minus) + sqrt(meas.PD_plus) );

% return to initilaized laser detuning value
if chl == 1
    AG4.DC(chl, prm.WP.lasers.pump_detuning);
end
if chl == 2
    AG4.DC(chl, prm.WP.lasers.probe_detuning);
end

%Return scope channels to Probe of 1:1
curr_inst.Probe(chn_minus, 1);
curr_inst.Probe(chn_plus, 1);


%Set initial Vscales
curr_inst.setVscale(chn_diff, Vscale_initial(1));
curr_inst.setVscale(chn_minus, Vscale_initial(2));
curr_inst.setVscale(chn_plus, Vscale_initial(3));

% vertical offsets
curr_inst.setVoffset(chn_diff, 0);
curr_inst.setVoffset(chn_minus, 0);
curr_inst.setVoffset(chn_plus, 0);
        
if flg.plot_data

    figure;
%     subplot(211)
    hold on;
    plot(meas.laser_det, meas.PD_abs, 'o-', 'markerSize', 4);
    plot(meas.laser_det, meas.PD_minus, 'v-', 'markerSize', 4);
    plot(meas.laser_det, meas.PD_plus, '^-', 'markerSize', 4);
    plot(meas.laser_det, meas.PD_diff, 'd-', 'markerSize', 4);
    legend('Abs.', 'PD -', 'PD +', 'PD diff.');
    xlabel('Laser detuning (volts)');
    ylabel('PD voltage (V)');
    title(['Rubidium absorption vs ' prm.which_laser ' laser detuning']);
    hold off;
%     subplot(212)
%     plot(meas.laser_det, meas.Rb_a85, '.', 'markerSize', 10)
%     title(['Rubidium 85 amplitude vs ' prm.which_laser ' laser detuning'])
%     xlabel('Laser detuning (volts)')
%     ylabel('Amplitude (volts)')
%     legend('Amplitude', 'Location', 'Best');

end

logger.log('=======  FID rubidium vs Laser wavelength experiment is done.  =======')
