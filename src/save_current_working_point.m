% ========================================================================
% This script reads the current working point from the lab instruments
% ========================================================================


%% TODO
% Add ESR locking LIA readout 19/09/21 - Roy
% Logging of Laser currents and cell temperature - CF


%% coils factors

%% Read current work point from instruments
% Fields
WP.fields.Bz.B0 = 80e-3; %DorB.Bz.apply('getCurrent');
WP.fields_control.Bz = 'Manual';

fields = {'Bx','By'};
for ind = 1:2
    field = fields{ind};
    [~,wvfm] = DorB.(field).apply('Output',[],true);
    WP.fields.(field).B0 = wvfm{1}.p(end);
    WP.fields.(field).Wfm = wvfm{1}.name;
    WP.fields.(field).params =  wvfm{1}.p;
    if ~strcmpi(wvfm{1}.name,'DC')
        [~,sync] = DorB.(field).apply('Syncio');
        WP.fields.(field).sync = sync;
    else
        WP.fields.(field).sync = 0;
    end
    WP.fields_control.(field) = 'Voltage';
end

fields = {'gradBx','gradBy','gradBz'};
for ind = 1:3
    field = fields{ind};    
    if isfield(DorB,field)
        WP.fields.(field) = DorB.(field).apply('getCurrent');
        WP.fields_control.(field) = 'Current';
    else
        WP.fields.(field) = 0;
        WP.fields_control.(field) = 'Manual';
    end
end
 

% Drives, LIAs and PIDs
drives = {'ESR','Xe129','Xe131'};
for ind=1:3
    drive = drives{ind};
    [stt,drive_params] = DorB.([drive 'Drive']).apply('Output',[],true);
    WP.(drive).drive_state = stt;
    WP.(drive).freq = drive_params{1}.p(1);
    WP.(drive).amp = drive_params{1}.p(2);
    [~,sync] = DorB.([drive 'Drive']).apply('Syncio');
    WP.(drive).sync = sync;

    if isfield(DorB,[drive '_LIA'])
        WP.(drive).LIAconfig = DorB.([drive '_LIA']).apply('readConfiguration');
        WP.(drive).LIA_ENBW = DorB.([drive '_LIA']).apply('ENBW');
    end

    [stt, modDev] = DorB.([drive 'Drive']).apply('FMmodulation');
    WP.(drive).mod_state = stt;
    WP.(drive).mod_dev = modDev;
    if isfield(DorB,[drive '_PID'])
        WP.(drive).PIDconfig = DorB.([drive '_PID']).apply('readConfiguration');
    end
end

fields = {'Bx','By'};
for ind = 1:2
    field = fields{ind};
    if isfield(DorB,[field '_sim'])
        [stt,wvfm] = DorB.([field '_sim']).apply('Output',[],true);
        phase = DorB.([field '_sim']).apply('getPhase');
        WP.Sim.(field).state = stt;      
        WP.Sim.(field).params =  num2cell([wvfm{1}.p(1:2) phase wvfm{1}.p(3)]);
    end
end
if isfield(DorB,'Sim_LIA')
    WP.Sim.LIAconfig = DorB.Sim_LIA.apply('readConfiguration');
end


% Trigger
[sst,wvfm] = DorB.Trigger.apply('Output',[],true);
WP.trigger.state = stt;
WP.trigger.Wfm = wvfm{1}.name;
WP.trigger.freq = wvfm{1}.p(1);
WP.trigger.amp = wvfm{1}.p(2);
WP.trigger.offset = wvfm{1}.p(3);


% Lasers
[~,wvfm] = DorB.LaserTempPump.apply('Output',[],true); WP.lasers.pump.temperature  = wvfm{1}.p(end);
[~,wvfm] = DorB.LaserTempProbe.apply('Output',[],true); WP.lasers.probe.temperature = wvfm{1}.p(end);
WP.lasers.pump.current  = 169e-3;      
WP.lasers.probe.current = 145.8e-3;
if isfield(DorB,'DeltaPump')
    WP.lasers.pump.detuning = DorB.DeltaPump.getDetuning;
end
if isfield(DorB,'DeltaProbe')
    WP.lasers.probe.detuning = DorB.DeltaProbe.getDetuning;
end


% Heater
% WP.T_cell = 108;%str2double(inputdlg('enter current temperature [C]:', 'temp', 1, {'108'}));

%%
prm.WP = WP;

logger.log('Reading current WP from instruments.')
clear WP
