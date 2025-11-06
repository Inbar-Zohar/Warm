%Settings:
measT = prm.cDAQ_meas_time_min; % min
ch_vec = prm.cDAQ_channels; % A vector of the channels we want to save
%Use cDAQ.getChannels() to get a list of available channels
Fs = prm.cDAQ_Fs;
prm.cDAQ_PAUSE_TIME_SEC = 5;

% experiment duration calculation
tot_time_sec = measT*60 + prm.cDAQ_PAUSE_TIME_SEC;
now_time = datetime;
end_time_hours = floor((tot_time_sec) / 60 / 60);
end_time_min = floor((tot_time_sec) / 60 - end_time_hours * 60);

if flg.is_part_of_loop
    fprintf('\n============================================================')
    fprintf(['\nStart Time:  ' datestr(datetime('now'))]);
    fprintf(['\nFinish Time: ' datestr(datetime('now') + duration([measT/60 + prm.cDAQ_PAUSE_TIME_SEC / 3600, 0, 0]))]);
    fprintf(['\nTotal Duration: ', num2str(end_time_hours),' hours ', num2str(end_time_min),' minutes. ']);
    fprintf('\n============================================================\n')
end
logger.log(['Total Duration: ', num2str(end_time_hours),' hours ', num2str(end_time_min),' minutes. '])

if Fs >51200
    Fs = 51200;
    warning('Fs was set to its maximal value of 51200 samp/sec');
end

%We assume the cDAQ object exists and is initialized!

%Clean all channels and set only the requested channels
% clean channels
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


cDAQ.setFs(Fs);

% clearvars meas;

% logger.log('Read cDAQ experiment.')
% if ~flg.is_part_of_loop
%     disp('[+] Read cDAQ experiment');
% end

%Get the channel information
chanInfo = cDAQ.getChannels();

%Set the measurements time
cDAQ.setDuration(measT*60); %Measurement time in seconds
if ~flg.is_part_of_loop
    logger.log('Reading DAQ...');
end

[t, v] = cDAQ.measureForeground; %Measure
% logger.log(('[+] done reading DAQ');
% pause(prm.cDAQ_PAUSE_TIME_SEC); %Pause for an extra time?

% tmpfname = 'cDAQTemp.mat';
% 

% 
% % 1) Set the cDAQ to measure for X min, and save data to a temporary file.
% lh = cDAQ.measureToFile(tmpfname, measT, false);
% 
% % 2) Wait X min
% pause(measT*60 + 5); %Extra 5 sec
% delete(lh); % Delete the listener handle
% cDAQ.Stop(); % And stop the cDAQ
% 
% % 3) Load the data from the file to the workspace
% load(tmpfname);
% % 4) Delete the file
% eval(['delete ' tmpfname ';']);
%Data is available in the workspace in the meas structure
meas=create_fig(meas); 
meas.cDAQ_t = t;
meas.cDAQ_v = v;
meas.cDAQ_chInfo = chanInfo;
% meas.lgd_cell = lgd_cell;

clearvars t v chanInfo;

% 5) Plot the data.

plot(meas.cDAQ_t,meas.cDAQ_v);
chanInfo = meas.cDAQ_chInfo;
lgd_cell = {chanInfo(:).Name};
legend(lgd_cell);
xlabel('Time [s]');
ylabel('Amplitude [V]');

meas.lgd_cell = lgd_cell;


%Finally, reset the cDAQ to Fs = 2048 samp/sec, and reconnect all channels
cDAQ.setFs(2048);
get_daq;
