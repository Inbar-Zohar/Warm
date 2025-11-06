% Nitzan Haviv 30.11.23

%Example code: in the same directory -> TEST.m
classdef ADP303XE
    
% The adderss parameter sould be:
% For tcpip - the ip address
% For usb   - Choose the X# from the list:
%    1#  Makat 2335811    ->         address = '0xF4EC::0x1430::SPD3XJFX7R1727'
%    2#  Makat 2335810    ->         address = '0xF4EC::0x1430::SPD3XJFX7R1728'

%For more information see:
%Team members/ Nitzan / connectivity / ADP303XE /Read me for USB connection

%Example: ADP_1 = ADP303XE('tcpip', '10.0.0.31'); 
%         ADP_2 = ADP303XE('usb', 2);
   
    properties
        Ins     % The connection to the PC
    end
    
    methods
        % Constructor
        function obj = ADP303XE(connection_type, address)
            if nargin > 0
                if ~exist('connection_type', 'var')
                    error('Please state the type of connection: tcpip / usb')
                end
                switch lower(connection_type)
                    case 'tcpip'
                        if isnumeric(address); address = num2str(address); end
                        	obj.Ins = visa('NI',strcat('TCPIP0::', address ,'::inst0::INSTR'));
                            fopen(obj.Ins);
                    case 'usb' % No need to specify the address
                        switch address
                            case 1
                                address = '0xF4EC::0x1430::SPD3XJFX7R1727';
                            case 2
                                address = '0xF4EC::0x1430::SPD3XJFX7R1728';
                            otherwise
                                Error('For USB connection, choose the device number (see comments)')
                        end
                        obj.Ins = visa('ni', ['USB0::' address '::INSTR']);
                        fopen(obj.Ins);
                end
            else
                error('Please specify the connection type. For tcpip specify the ip address. ');
            end
        end
        
        
        function data = IDN(obj) %Indentify
        command = '*IDN?' ;
        fprintf(obj.Ins, command);
        data = fscanf(obj.Ins,'%c');
        disp(data);
        end
        
        %% set functions: 
 
        function OutputON(obj, channel) 
            fprintf(obj.Ins, ['OUTPut CH' num2str(channel) ',ON']);
        end
        
        function OutputOFF(obj, channel)
            fprintf(obj.Ins, ['OUTPut CH' num2str(channel) ',OFF']);
        end
        
        function OutputOFF_All(obj)
            fprintf(obj.Ins, 'OUTPut CH1,OFF');
            fprintf(obj.Ins, 'OUTPut CH2,OFF');
            fprintf(obj.Ins, 'OUTPut CH3,OFF');
        end
       
        function setParallelMode(obj)
            command = 'OUTPut:TRACK 2';
            fprintf(obj.Ins, command);
        end
        
        function setSeriesMode(obj)
            command = 'OUTPut:TRACK 1';
            fprintf(obj.Ins, command);
        end
        
        function setIndependentMode(obj)
            command = 'OUTPut:TRACK 0';
            fprintf(obj.Ins, command);
        end
        
        function setCurrent(obj, channel, current)
            fprintf(obj.Ins, ['CH' num2str(channel) ':CURRent ' num2str(current)]);
        end
        
        function setVoltage(obj, channel, voltage)
            fprintf(obj.Ins, ['CH' num2str(channel) ':VOLTage ' num2str(voltage)]);
        end
         
        
        %% Get functions: 
                
        function response = getCurrent(obj, channel)
        command = 'MEASure:CURRent? CH' ;
        fprintf(obj.Ins, [command num2str(channel)]);
        response = str2double(fscanf(obj.Ins));
        end
        
        function response = getSetCurrent(obj, channel) %Returns the defined current for the channel
        fprintf(obj.Ins, ['CH' num2str(channel) ':CURRent?']);
        response = str2double(fscanf(obj.Ins));
        end
           
        function response = getVoltage(obj, channel)
        command = 'MEASure:VOLTage? CH' ;
        fprintf(obj.Ins, [command num2str(channel)]);
        response = str2double(fscanf(obj.Ins));
        end
        
        function response = getSetVoltage(obj, channel) %Returns the defined voltage for the channel
        fprintf(obj.Ins, ['CH' num2str(channel) ':VOLTage?']);
        response = str2double(fscanf(obj.Ins));
        end
        
        function response = getPower(obj, channel)
        command = 'MEASure:POWEr? CH' ;
        fprintf(obj.Ins, [command num2str(channel)]);
        response = str2double(fscanf(obj.Ins));
        end
        
        function response = getIP_address(obj)
        fprintf(obj.Ins, 'IPaddr?');
        response = fscanf(obj.Ins);
        disp(response);
        end
    

    
    
    
    end
end
