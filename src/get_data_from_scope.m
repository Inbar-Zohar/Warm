logger.log('Read Scope experiment.')
% if ~flg.is_part_of_loop
%     disp('[+] Read Scope experiment');
% end

if ( ischar(prm.scope_obj) && exist(prm.scope_obj,'var') )

scope_obj = evalin('base', prm.scope_obj);
    
channels = prm.channels_scope;
% gain = 10 / 20e-3;
% under_sample = 100000;
gains = prm.gains_scope;
% Sampling rate depends on the number of active channels. The scope has
% maximal memory of 62500 samples (when all channels are active without high res). Sampling
% rate in this case is Fs = 62500/Tscale/10;
% Fs = prm.Fs_scope; % Sp/sec
% Tscale = prm.Tscale_scope; % sec/div
%Get the scope time scale (sec/div)
Tscale = scope_obj.getTscale;

% Tdelay = 0;
% scope_obj.setTscale(Tscale);
% scope_obj.setTdelay(Tdelay)
% scope_obj.setTref('CENTer');
% scope_obj.TrigSource('EXT');
scope_obj.Single;
% pause(5 * Tscale);
% scope_obj.highres;

% start measurement
scope_obj.ForceTrig;
scope_obj.readyToRead;


% read signal from scope
[t_, v_] = scope_obj.Read(channels);
% t = t_(1:round(length(t_) / under_sample):end);
% v = v_(:, 1:round(length(v_) / under_sample):end);
t = t_;
v = v_;

% release scope
% scope_obj.setTref('LEFT'); 
% scope_obj.TrigSource('LINE'); 
scope_obj.Run;

% set local variables
meas.scope_experiment.time = t;
for i = 1:length(channels)
    if channels(i) == 1
        meas.scope_experiment.cha1 = v(i,:) / gains(i);
    end
    if channels(i) == 2
        meas.scope_experiment.cha2 = v(i,:) / gains(i);
    end
    if channels(i) == 3
        meas.scope_experiment.cha3 = v(i,:) / gains(i);
    end
    if channels(i) == 4
        meas.scope_experiment.cha4 = v(i,:) / gains(i);
    end
end

figure;
plot(t, v(1,:));
hold on;
lgd_cell = {};
lgd_cell{1} = num2str(channels(1));
for i = 2:length(channels)
    plot(t, v(i, :));
    lgd_cell{i} = num2str(channels(i));
end
grid on
legend(lgd_cell);
xlabel('Time [s]');
ylabel('Amplitude [V]');

%Tidy up
clearvars scope_obj t_ t v_ v;

else
    %Requested scope object doesn't exist
    error('Read Scope experiment faild: Requested scope object does not exist');
end

% 1/(t_(2)-t_(1))