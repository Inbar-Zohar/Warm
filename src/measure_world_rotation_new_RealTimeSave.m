% ========================================================================
%   This script implements a World Rotation measurement and lauches the
%   analysis script.
% ========================================================================

%% Experiment duration calculation
tot_time_sec = floor(60 * prm.exp_time_min) + prm.cdaq_pause_time_sec;
duration_hours = floor((tot_time_sec) / 60 / 60);
duration_min = floor((tot_time_sec) / 60 - duration_hours * 60);
logger.log(['Total Duration: ', num2str(duration_hours),' hours ', num2str(duration_min),' minutes'])
logger.log(['Finish Time: ', datestr(datetime('now') + duration([duration_hours, duration_min, 0]))])

%% Verify PID locks
if flg.verify_pid_locks
    %%
    Xenons = {'Xe129','Xe131'};
    for ind=1:2
        Xe = Xenons{ind};
        LIA = [Xe '_LIA'];
        PID = [Xe '_PID'];
        if isfield(DorB,LIA) && isfield(DorB,PID)
            % if ~DorB.(PID).apply('getOutputMode')
                DorB.(LIA).apply('AutoPhase')
                pause(0.5)
                DorB.(PID).apply('setOutputModePID')
            % end
        else
            warning('Cannot verify PID lock, devices not defined.')
        end
    end
end

%% read DAQ
pause(prm.cdaq_pause_time_sec);
if prm.Fs > 51200
    prm.Fs = 51200;
    warning('Fs was set to its maximal value of 51200 samp/sec');
end
if prm.Fs < 2048
    prm.Fs = 2048;
    warning('Fs was set to its minimal value of 2048 samp/sec');
end
cDAQ.setFs(prm.Fs);
meas.chanInfo = cDAQ.getChannels();
meas.chan_names = {meas.chanInfo.Name};
cdaq_filename = 'cdaq_meas.mat';

cDAQ.setDuration(tot_time_sec);
logger.log('Reading DAQ and saving files...');
lh = cDAQ.measureToFile(cdaq_filename,prm.save_time_min,false);
clear vpr
clear tpr
for ind=[1,2]
    logger.log(sprintf('saving priod %i',ind))
    
    tries = 0; success = 0;
    while tries < tot_time_sec/10 && ~success
        try
            cdaq_data = load(cdaq_filename);
            success = 1;
        end
        if ~success
            tries = tries + 1;
            logger.log(sprintf('still waiting %i',tries))
            pause(11);
        end
    end
    if ~success; error('Didn''t get cdaq file'); end
    
    vpr_temp = cdaq_data.v;
    tpr_temp = cdaq_data.t;
    if ind==1
        tpr = tpr_temp';
        vpr = vpr_temp';
    else
        tpr = cat(2,tpr,tpr_temp');
        vpr = cat(2,vpr,vpr_temp');
    end
    %% Apply requested smoothing
    if prm.data_smooth_param > 1
        vpr_smoothed = decimate_array(vpr,prm.data_smooth_param);
        tpr_smoothed = tpr(floor(prm.data_smooth_param/2):prm.data_smooth_param:end);
        if length(tpr_smoothed) > size(vpr_smoothed,2); tpr_smoothed = tpr_smoothed(1:size(vpr_smoothed,2)); end %PROBABLY UNNECESSARY
        logger.log(sprintf('Smoothed data by %i samples',prm.data_smooth_param));
    else
        tpr_smoothed = tpr;
        vpr_smoothed = vpr;
    end
    
    %% Check data size and apply extra smoothing if relevant
    meas_info = whos('vpr_smoothed');
    size_GB = meas_info.bytes / (2^3)^10;
    if size_GB > prm.max_full_variable_size_GB
        if flg.allow_variables_above_max
            warning('Full variable size is %.1f', size_GB);
        else
            extra_smooth_param = ceil(size_GB / prm.max_full_variable_size_GB);
            warning('Extra smoothing by %i samples to reduce full variable size below 2GB, might lose information', extra_smooth_param)
            vpr_smoothed = decimate_array(vpr_smoothed,extra_smooth_param);
            tpr_smoothed = tpr_smoothed((extra_smooth_param/2):extra_smooth_param:end);
            if length(tpr_smoothed) > size(vpr_smoothed,2); tpr_smoothed = tpr_smoothed(1:size(vpr_smoothed,2)); end %PROBABLY UNNECESSARY
            logger.log(sprintf('Smoothed data by %i samples',extra_smooth_param));
        end
    end
    
    %% Verify the PID channels are actually saved 
    ch_names = {meas.chanInfo.Name};
    PID131_ch = find(strcmp('PID-131 out', {meas.chanInfo.Name}),1);
    if ~isempty(PID131_ch) && ~ismember(PID131_ch,prm.ch_to_save)
        prm.ch_to_save = [prm.ch_to_save PID131_ch];
        warning('It seems you didn''t want to save the 131 PID output channel. Well, nobody''s asking you!')
    end
    PID129_ch = find(strcmp('PID-129 out', {meas.chanInfo.Name}),1);
    if ~isempty(PID129_ch) && ~ismember(PID129_ch,prm.ch_to_save)
        prm.ch_to_save = [prm.ch_to_save PID129_ch];
        warning('It seems you didn''t want to save the 129 PID output channel. Well, nobody''s asking you!')
    end
    
    %% Extract data to fields
    meas.time_s = tpr_smoothed;
    for ind=prm.ch_to_save
        this_name = strrep(cDAQ.Ins.Channels(ind).Name,' ','_');
        this_name = strrep(this_name, '-', '_');
        if isempty(this_name); this_name = 'Empty'; end
        meas.(this_name) = vpr_smoothed(ind,:);    
    end
    
    %% Analyze, save and plot
    if flg.keep_raw_data
        raw_meas_info = whos('vpr');
        size_GB = raw_meas_info.bytes / (2^3)^10;
        if size_GB > prm.max_full_variable_size_GB
            warning('Saving raw data even though the size is %.1f GB, you crazy person',size_GB);
        else
            warning('Saving raw data. Size isn''t too large but this should be done only for debugging');        
        end
        meas.t_raw = tpr;
        meas.v_raw = vpr;
    end
    
    logger.log('Extracting World Rotation')
    [Omega_deg_hr,Bz,Omega_129_deg_hr,Omega_131_deg_hr] = ...
        get_Omega_and_Bz_from_signals_new(meas,prm);
    meas.Omega_deg_hr = Omega_deg_hr;
    meas.Bz = Bz;
    if flg.plot_individual_xenons
        meas.Omega_129_deg_hr = Omega_129_deg_hr;
        meas.Omega_131_deg_hr = Omega_131_deg_hr;
    end
    meas.fs = 1/diff(meas.time_s(1:2));
    
    if flg.plot_measured_signals
        try
            plot_world_rotation_measured_signals;
        end
    end
    
    % Get Noise Spectrum
    if flg.analyze_rotation_signal
        logger.log('Analyzing World Rotation')
        if isfield(DorB,'Xe129_LIA') && isfield(DorB,'Xe131_LIA')
            prm.estimated_bw = min([prm.estimated_bw DorB.Xe129_LIA.apply('ENBW') DorB.Xe131_LIA.apply('ENBW')]);
        end
        analyze_world_rotation;    
    end
    delete(cdaq_filename);
end
% Plot data
% if flg.plot_data
%     %%
%     logger.log('Plotting raw data')
%     meas = create_fig(meas,1010);
%     % plot_single_experiment(meas, prm);
%     plot(tpr_smoothed,vpr_smoothed);
%     lgd_cell = {meas.chanInfo(:).Name};
%     legend(lgd_cell);
%     xlabel('Time [s]');
%     ylabel('Amplitude [V]');
% end

%% Reset the cDAQ
% clear tpr vpr tpr_smoothed vpr_smoothed
cDAQ.clearListener(lh);
get_daq;
