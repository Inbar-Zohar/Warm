% ========================================================================
% This script loads a working point to the lab instruments from file
% ========================================================================


%% TODO

%% Open folder from which to load the Working Point
address = save_path;
[fname, address] = uigetfile(address,'Select WorkingPoint to load to system');
temp = load([address fname], 'prm');
prm.WP = temp.prm.WP;
WP = prm.WP;

%% Write current work point to instruments
if ~isfield(WP,'fields_control')
    warning('Old file, things may fail.')
end

%% Fields

% Bx, By
fields = {'Bx','By'};
for ind = 1:2
    field = fields{ind};
    if strcmp(WP.fields_control.(field),'Voltage')
        if isfield(WP.fields,field)
            is_set = 0;
            if strcmpi(WP.fields.(field).Wfm,'SIN')
                params = (WP.fields.(field).params);
                params = num2cell([params(1:2) 0 params(3)]);
                DorB.(field).apply('Sin',params{:});
                if WP.fields.(field).sync; DorB.(field).apply('Syncio',1); end
                is_set = 1;
            else
                DorB.(field).apply('DC',WP.fields.(field).B0);
                is_set = 1;
            end
            DorB.(field).apply('OutputON')
            disp([field ' set']);
        elseif ~isfield(WP.fields,field) && isfield(WP.fields,[field '0'])
            DorB.(field).apply('DC',WP.fields.([field '0']));
            DorB.(field).apply('OutputON')
            disp([field '0 set']);
        end
    end
end

% Bz
if isfield(DorB,'Bz') && isfield(WP.fields, 'Bz0') 
    if strcmp(WP.fields_control.Bz,'Current')
        DorB.Bz.apply('setCurrent', WP.fields.Bz0);  DorB.Bz.apply('OutputON');
    end
    disp('Bz set');        
end

fields = {'gradBx','gradBy','gradBz'};
for ind = 1:3
    field = fields{ind};
    if isfield(DorB, field) && isfield(WP.fields, field)
        if strcmp(WP.fields_control.(field),'Current')
            DorB.(field).apply('OutputON');
            DorB.(field).apply('setCurrent', WP.fields.(field));
        end
        disp([field ' set']);
    end
end


%% Drives
drives = {'ESR','Xe129','Xe131'};
for ind=1:3
    drive = drives{ind};
    drive_dev = ([drive 'Drive']);
    DorB.(drive_dev).apply('Sin',WP.(drive).freq,WP.(drive).amp);
    DorB.(drive_dev).apply('Output',WP.(drive).drive_state);
    if WP.(drive).sync; DorB.(drive_dev).apply('Syncio',1); end

    DorB.([drive 'Drive']).apply('FMmodulation',WP.(drive).mod_state,WP.(drive).mod_dev);

    if isfield(DorB,[drive '_LIA']) && isfield(WP.(drive),'LIAconfig')
        DorB.([drive '_LIA']).apply('applyConfiguration',WP.(drive).LIAconfig);
    end

    if isfield(DorB,[drive '_PID']) && isfield(WP.(drive),'PIDconfig')
        DorB.([drive '_PID']).apply('applyConfiguration',WP.(drive).PIDconfig);
    end

    disp([drive ' set']);
end

if isfield(WP,'Sim')
    fields = {'Bx','By'};
    for ind = 1:2
        field = fields{ind};
        if isfield(DorB,[field '_sim']) && isfield(WP.Sim,field)
            DorB.([field '_sim']).apply('Sin',WP.Sim.(field).params{:});
            if WP.Sim.(field).state 
                DorB.([field '_sim']).apply('OutputON');
            else
                DorB.([field '_sim']).apply('OutputOFF');
            end                
        end
    end
    if isfield(DorB,'Sim_LIA') &&  isfield(WP.Sim,'LIAconfig')
        DorB.Sim_LIA.apply('applyConfiguration',WP.Sim.LIAconfig);
    end
    disp('Simulation set');        
end

%% Trigger
switch WP.trigger.Wfm
    case 'SQU'
        funcname = 'Square';
    case 'SIN'
        funcname = 'Sin';
end
params = num2cell([WP.trigger.freq WP.trigger.amp WP.trigger.offset]);
DorB.Trigger.apply(funcname,params{:});
if WP.trigger.state
    DorB.Trigger.apply('OutputON');
else
    DorB.Trigger.apply('OutputOFF');
end

%% load the laser detuning only if flag is 1
if flg.load_laser_params
    % pump detuning
    DorB.LaserTempPump.apply('DC',lasers.pump.temperature); DorB.LaserTempPump.apply('OutputON')
    if isfield(DorB,'DeltaPump') && isfield(WP.lasers.pump,'detuning')
        DorB.DeltaPump.setDetuning(WP.lasers.pump.detuning)
    end
    disp('Pump detuning set');
        
    % probe detuning
    DorB.LaserTempProbe.apply('DC',lasers.probe.temperature); DorB.LaserTempProbe.apply('OutputON')
    if isfield(DorB,'DeltaProbe') && isfield(WP.lasers.probe,'detuning')
        DorB.DeltaProbe.setDetuning(WP.lasers.probe.detuning)
    end
    disp('Probe detuning set');
end
%%

logger.log(['WP loaded from file --- ' fullfile(address,fname)])
clear WP
