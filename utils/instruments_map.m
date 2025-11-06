function [inst_id, daq_id] = instruments_map()
%instruments_map is a map between lab instruments to physical objects
% measured.
    
    % Agilent address
    inst_id.AG1.serial = '0x2C07::MY58001471';
    inst_id.AG1.connection_type = 'usb';
    inst_id.AG2.serial = '192.168.1.41';
    inst_id.AG2.connection_type = 'tcpip';
    inst_id.AG3.serial = '0x2C07::MY58001472';
    inst_id.AG3.connection_type = 'usb';
    inst_id.AG4.serial = '192.168.1.42';
    inst_id.AG4.connection_type = 'tcpip';
    inst_id.AG6.serial = '192.168.1.43';
    inst_id.AG6.connection_type = 'tcpip';
    
    % Coils address
%     inst_id.BB.Address = '192.168.1.44';
%     inst_id.BK.com = 1;
    inst_id.CS1.com = 3;
    
    % ESR LIA address
    inst_id.SR_ESR.address = 30;
    inst_id.SR_ESR.GPIB_address = 7;
    inst_id.SR_129.address = 29;
    inst_id.SR_129.GPIB_address = 7;
    inst_id.SR_131.address = 28;
    inst_id.SR_131.GPIB_address = 7;
    
    % PIDs address
    inst_id.PIDs.com = 7;
    
    % DAQ address
    daq_id.probe_temperature.address = [2, 0];
    daq_id.probe_temperature.name = 'probe_temperature';
    daq_id.PID_131_out.address = [2, 1];
    daq_id.PID_131_out.name = 'Detuning-131 correction - PID-131 output';
    daq_id.PID_129_error.address = [1, 0];
    daq_id.PID_129_error.name = 'Xe-129 theta - PID-129 error';
    daq_id.PID_131_error.address = [1, 1];
    daq_id.PID_131_error.name = 'Xe-131 theta - PID-131 error';
    daq_id.Xe_129_R.address = [2, 2];
    daq_id.Xe_129_R.name = 'Xe-129 R';
    daq_id.Xe_131_R.address = [2, 3];
    daq_id.Xe_131_R.name = 'Xe-131 R';
    daq_id.Delta_Bz.address = [1, 2];
    daq_id.Delta_Bz.name = 'Delta Bz';
    daq_id.pump_temperature.address = [1, 3];
    daq_id.pump_temperature.name = 'pump_temperature';
end

