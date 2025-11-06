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
cDAQ.setDuration(tot_time_sec);
logger.log('Reading DAQ...');
[tpr, vpr] = cDAQ.measureForeground;
tpr = tpr';
vpr = vpr';
logger.log('Done reading DAQ');

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

% Plot data
if flg.plot_data
    %%
    logger.log('Plotting raw data')
    meas = create_fig(meas,1010);
    % plot_single_experiment(meas, prm);
    plot(tpr_smoothed,vpr_smoothed);
    lgd_cell = {meas.chanInfo(:).Name};
    legend(lgd_cell);
    xlabel('Time [s]');
    ylabel('Amplitude [V]');
end

%% Reset the cDAQ
clear tpr vpr tpr_smoothed vpr_smoothed
get_daq;
