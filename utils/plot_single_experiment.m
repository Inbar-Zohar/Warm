function plot_single_experiment(daq_data, exp_data)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    
    sens_V_Theta = 10/180; %V/deg
    color129 = [0.1 0 0.95];
    color131 = [0.8 0.2 0];
    
    subplot(3, 2, 2)
    hold on
    if isfield(daq_data,'Bz_V')
        plot(daq_data.time_s, daq_data.Bz_V, 'color', color129, 'displayname', ['dV Xe129'])
    end
    if isfield(daq_data,'Det131_PIDb_V')
        plot(daq_data.time_s, daq_data.Det131_PIDb_V, 'color', color131, 'displayname', ['dV Xe131'])
    end
    ylabel('[V]')
    xlabel('time [sec]')
    title('PID loop output')
    legend show
    subplot(3, 2, 4)
    hold on
    if isfield(daq_data,'Xe129_R_V')
        plot(daq_data.time_s, daq_data.Xe129_R_V, 'color', color129, 'displayname', ['Xe129 R'])
    end
    if isfield(daq_data,'Xe131_R_V')
        plot(daq_data.time_s, daq_data.Xe131_R_V, 'color', color131, 'displayname', ['Xe131 R'])
    end
    ylabel('[V]')
    xlabel('time [sec]')
    legend show
    subplot(3, 2, 6)
    hold on
    if isfield(daq_data,'Xe129_theta_V')
        plot(daq_data.time_s, daq_data.Xe129_theta_V / sens_V_Theta * exp_data.DAQ_ATT_129, 'color', color129, 'displayname', ['Xe129 \theta'])
    end
    if isfield(daq_data,'Xe131_theta_V')
        plot(daq_data.time_s, daq_data.Xe131_theta_V / sens_V_Theta * exp_data.DAQ_ATT_131, 'color', color131, 'displayname', ['Xe131 \theta'])
    end
    ylabel('[Deg]')
    xlabel('time [sec]')
    legend show 
    
    subplot(3, 2, [1, 3, 5])
    hold on
%     plot(daq_data.time_s, daq_data.ProbeT_V, 'displayname', ['ProbeT'])
    if isfield(daq_data,'ESR_LIA_out')
        plot(daq_data.time_s, daq_data.ESR_LIA_out, 'displayname', ['ESR LIA out V_{RMS}'])
    end    
if isfield(daq_data,'Det131_PIDb_V')
        plot(daq_data.time_s, daq_data.Det131_PIDb_V, 'displayname', ['Det131 PIDb'])
    end
    if isfield(daq_data,'Xe129_theta_V')
        plot(daq_data.time_s, daq_data.Xe129_theta_V, 'displayname', ['Xe129 \theta'])
    end
    if isfield(daq_data,'Xe131_theta_V')
        plot(daq_data.time_s, daq_data.Xe131_theta_V, 'displayname', ['Xe131 \theta'])
    end
    if isfield(daq_data,'Xe129_R_V')
        plot(daq_data.time_s, daq_data.Xe129_R_V, 'displayname', ['Xe129 R'])
    end
    if isfield(daq_data,'Xe131_R_V')
        plot(daq_data.time_s, daq_data.Xe131_R_V, 'displayname', ['Xe131 R'])
    end
    if isfield(daq_data,'Bz_V')
        plot(daq_data.time_s, daq_data.Bz_V, 'displayname', ['Bz'])
    end
%     plot(daq_data.time_s, daq_data.PumpT_V, 'displayname', ['PumpT'])
    
    title_str = '';
    if isfield(daq_data,'DateTime_var')
        title_str = [title_str char(daq_data.DateTime_var) '\n'];
    end
    title_str = [title_str exp_data.name];

    title(sprintf(title_str))
    legend show
    
end

