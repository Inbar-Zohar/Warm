% ========================================================================
%   This script implementing a World Rotation Measurement and computing the
%   power spectral density (PSD) from the measured data. It assumes that
%   the system is in an optimal working point.
% ========================================================================

% Initial setup
prm.to_plot = prm.points_to_plot;
prm.readtimeM = prm.exp_time_min;
prm.readtimeH = prm.readtimeM / 60;
prm.readtimeS = floor(60 * prm.readtimeM);
prm.name = [num2str(prm.readtimeM) ' mins'];

% experiment duration calculation
tot_time_sec = prm.readtimeS + prm.PAUSE_TIME_SEC;
now_time = datetime;
end_time_hours = floor((tot_time_sec) / 60 / 60);
end_time_min = floor((tot_time_sec) / 60 - end_time_hours * 60);
if ~flg.is_part_of_loop
    fprintf('\n============================================================')
    fprintf(['\nStart Time:  ' datestr(datetime('now'))]);
    fprintf(['\nFinish Time: ' datestr(datetime('now') + duration([prm.readtimeH + prm.PAUSE_TIME_SEC / 3600, 0, 0]))]);
    fprintf(['\nTotal Duration: ', num2str(end_time_hours),' hours ', num2str(end_time_min),' minutes. ']);
    fprintf('\n============================================================\n')
end
logger.log(['Total Duration: ', num2str(end_time_hours),' hours ', num2str(end_time_min),' minutes. '])
% logger.log('World rotation measurement initialized')
% if ~flg.is_part_of_loop
%     disp('[+] Starting a World Rotation Measurement');
% end

if prm.Fs >51200
    prm.Fs = 51200;
    warning('Fs was set to its maximal value of 51200 samp/sec');
end

%--- cDAQ -----------------------------
%We only need the "PID-131 out", and "PID-129 out" channels, so keep only
%them. Run cDAQ.getChannels() to see the list of available channels.

%Add additional channels in here if needed:
% prm.cDAQ_channels = [3 4]; %Only "PID-131 out", and "PID-129 out"
prm.cDAQ_channels = [1 2 3 4 5 6 7 8 9 10 11 12]; %Everything
if ~sum(prm.cDAQ_channels == 1)
    cDAQ.removeCh(1, 0);
end
if ~sum(prm.cDAQ_channels == 2)
    cDAQ.removeCh(1, 1);
end
if ~sum(prm.cDAQ_channels == 3)
    cDAQ.removeCh(1, 2);
end
if ~sum(prm.cDAQ_channels == 4)
    cDAQ.removeCh(1, 3);
end
if ~sum(prm.cDAQ_channels == 5)
    cDAQ.removeCh(2, 0);
end
if ~sum(prm.cDAQ_channels == 6)
    cDAQ.removeCh(2, 1);
end
if ~sum(prm.cDAQ_channels == 7)
    cDAQ.removeCh(2, 2);
end
if ~sum(prm.cDAQ_channels == 8)
    cDAQ.removeCh(2, 3);
end
if ~sum(prm.cDAQ_channels == 9)
    cDAQ.removeCh(3, 0);
end
if ~sum(prm.cDAQ_channels == 10)
    cDAQ.removeCh(3, 1);
end
if ~sum(prm.cDAQ_channels == 11)
    cDAQ.removeCh(3, 2);
end
if ~sum(prm.cDAQ_channels == 12)
    cDAQ.removeCh(3, 3);
end
%--------------------------------------


% get values from instruments
[~, sys_data.DATAT129] = AG1.Output(1, [], true); sys_data.DATAT129 = sys_data.DATAT129{1}.p(1:2);
[~, sys_data.DATAT131] = AG2.Output(1, [], true); sys_data.DATAT131 = sys_data.DATAT131{1}.p(1:2);
[~, sys_data.DATATESR] = AG3.Output(1, [], true); sys_data.DATATESR = sys_data.DATATESR{1}.p(1:2);
[~, sys_data.DATATVBx] = AG1.Output(2, [], true); sys_data.DATATVBx = sys_data.DATATVBx{1}.p(end);
[~, sys_data.DATATVBy] = AG2.Output(2, [], true); sys_data.DATATVBy = sys_data.DATATVBy{1}.p(end);
[~, sys_data.detPump ] = AG4.Output(1, [], true); sys_data.detPump  = sys_data.detPump{1}.p(end);
[~, sys_data.detProbe] = AG4.Output(2, [], true); sys_data.detProbe = sys_data.detProbe{1}.p(end);

% Bias fields
prm.Bx = sys_data.DATATVBx;  % [V]
prm.By = sys_data.DATATVBy;  % [V]

% ESR data
prm.fresESR = sys_data.DATATESR(1);  %[Hz]
prm.AmpESR = [sys_data.DATATESR(2)]; %[Vpp]

% NMR Xe129
prm.fres129 = sys_data.DATAT129(1);   %[Hz]
prm.Amp129 = [sys_data.DATAT129(2)];  %[Vpp]

% NMR Xe131
prm.fres131 = sys_data.DATAT131(1);   %[Hz]
prm.Amp131 = [sys_data.DATAT131(2)];  %[Vpp]

% Modulation
prm.df131 = str2double(query(AG2.Ins,'SOURce1:FM:DEV?')); %[Hz]
if prm.american_lock
    prm.df129 = str2double(query(AG1.Ins,'SOURce1:FM:DEV?')); %[Hz]
else
    prm.df129 = prm.df131 ; %[Hz]
end

% Xenons = {'Xe129','Xe131'};
% for ind=1:2
%     Xe = Xenons{ind};
%     LIA = [Xe '_LIA'];
%     PID = [Xe '_PID'];
%     if isfield(DorB,LIA) && isfield(DorB,PID)
%         if ~DorB.(PID).apply('getOutputMode')
%             DorB.(LIA).apply('AutoPhase')
%             pause(0.5)
%             DorB.(PID).apply('setOutputModePID')
%         end
%     end
% end

% read DAQ
pause(prm.PAUSE_TIME_SEC);
cDAQ.setFs(prm.Fs);
meas.chanInfo = cDAQ.getChannels();
cDAQ.setDuration(prm.readtimeS + prm.PAUSE_TIME_SEC);
logger.log('Reading DAQ...');
[tpr, vpr] = cDAQ.measureForeground;
logger.log('Done reading DAQ');
% if the data is too large, reduce to single precision and take 1/10th of
% the data
% if numel(vpr)*8>7.078051776e+09/10*2 % larger than 2 GB
% 
%     vpr_tmp = vpr(1:10:end,:);
%     vpr_sngl = single(vpr_tmp);
%     meas.tp_wrld_rot = tpr(1:10:end)';
%     meas.vp_wrld_rot = vpr_sngl';
%     clear vpr_tmp vpr_sngl
% else
meas.tp_wrld_rot = tpr';
meas.vp_wrld_rot = vpr';
% end
clear sys_data

%% unpack and smooth data
if 0
    for I=1:numel(prm.cDAQ_channels)
        %     this_data = smooth(meas.vp_wrld_rot(I, :), prm.data_smooth_param);
        this_data = meas.vp_wrld_rot(I, :);
        warning('removed smoothing!!!!!!!!!!!!!!!')
        switch prm.cDAQ_channels(I)
            case 1
                %ESR LIA out - RAW
                ESR_LIA_out   = this_data;
                meas.ESR_LIA_out = ESR_LIA_out(1:prm.data_smooth_param:end)';
            case 2
                %ESR LIA out - HP filtered
            case 3
                %PID-131 out
                Det131_PIDb_V   = this_data;
                meas.Det131_PIDb_V = Det131_PIDb_V(1:prm.data_smooth_param:end)';
            case 4
                %PID-129 out
                Bz_V            = this_data;
                meas.Bz_V = Bz_V(1:prm.data_smooth_param:end)';
            case 5
                %Xe-129 theta - PID-129 error
                Xe129_theta_V   = this_data;
                meas.Xe129_theta_V = Xe129_theta_V(1:prm.data_smooth_param:end)';
            case 6
                %Xe-131 theta - PID-131 error
                Xe131_theta_V   = this_data;
                meas.Xe131_theta_V = Xe131_theta_V(1:prm.data_smooth_param:end)';
            case 7
                %Xe-129 R
                Xe129_R_V       = this_data;
                meas.Xe129_R_V = Xe129_R_V(1:prm.data_smooth_param:end)';
            case 8
                %Xe-131 R
                Xe131_R_V       = this_data;
                meas.Xe131_R_V = Xe131_R_V(1:prm.data_smooth_param:end)';
            case 11
                %Xe-131 R
                Probe_after_cell       = this_data;
                meas.Probe_after_cell = Probe_after_cell(1:prm.data_smooth_param:end)';
                
        end
    end
    clear ProbeT_V Det131_PIDb_V Xe129_theta_V Xe131_theta_V Xe129_R_V Xe131_R_V Bz_V PumpT_V;

else
    %%
    for I=prm.ch_to_save
        this_name = strrep(cDAQ.Ins.Channels(I).Name,' ','_');
        this_name = strrep(this_name, '-', '_');
        this_data = meas.vp_wrld_rot(I, :);
        
        % eval([this_name ' = this_data;']);
        % eval(['meas.' this_name ' = smooth(this_data(1:prm.data_smooth_param:end));']);       
        eval(['meas.' this_name ' = decimate(this_data,prm.data_smooth_param);']);

        switch I 
            case 3
                eval(['meas.Det131_PIDb_V = meas.' this_name ';']);
            case 4
                eval(['meas.Bz_V = meas.' this_name ';']);
        end
    end
end
meas.time_s = meas.tp_wrld_rot(floor(prm.data_smooth_param/2):prm.data_smooth_param:end);

if isfield(meas,'tp_wrld_rot'); meas=rmfield(meas,'tp_wrld_rot'); end
if isfield(meas,'vp_wrld_rot'); meas=rmfield(meas,'vp_wrld_rot'); end

% ProbeT_V        = smooth(meas.vp_wrld_rot(1, :), prm.data_smooth_param);
% PumpT_V         = smooth(meas.vp_wrld_rot(8, :), prm.data_smooth_param);

% meas.ProbeT_V = ProbeT_V(1:prm.data_smooth_param:end)';
% meas.PumpT_V = PumpT_V(1:prm.data_smooth_param:end)';
%%


% Get Noise Spectrum
if flg.get_noise_spectrum
%     logger.log('Compute noise spectrum from data')
%     %[f,N,T,S,tau] = get_noise_spectrum_from_PID_output(meas, prm);
%     meas = create_fig(meas, 113);
%     clf;
% 
%     [f,N,T,S,tau, noise_mdl_prm] = get_noise_spectrum_with_fit(meas, prm, prm.estimated_bw);
%     %     [f,N,T,S,tau, noise_mdl_prm] = get_noise_spectrum_with_fit_kostya(meas, prm);
%     meas.noise_mdl_prm = noise_mdl_prm;
%     
%     meas.dBz = S(1, :); 
%     meas.Omega_deg_hr = S(2,:);
    logger.log('Analyzing World Rotation')
    if isfield(DorB,'Xe129_LIA') && isfield(DorB,'Xe131_LIA')
        prm.estimated_bw = min([prm.estimated_bw DorB.Xe129_LIA.apply('ENBW') DorB.Xe131_LIA.apply('ENBW')]);
    end
    analyze_world_rotation;
end

% Plot data
if flg.plot_data
    logger.log('Plot data')
    meas = create_fig(meas,1010);
    % plot_single_experiment(meas, prm);
    plot(meas.tp_wrld_rot,meas.vp_wrld_rot);
    lgd_cell = {meas.chanInfo(:).Name};
    legend(lgd_cell);
    xlabel('Time [s]');
    ylabel('Amplitude [V]');
end

%Reset the cDAQ
get_daq;

% logger.log('World Rotation measurement done.')
