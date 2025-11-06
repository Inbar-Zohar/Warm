% Nitzan Haviv 02.08.23

%Example code: in the same directory -> TEST.m
classdef CS580
    
    % This device, uses RS-232 interface (to connect the PC).

    properties
        Ins % The COM port used to connect the device
    end
    
    methods
        % Constructor
        function obj = CS580(com)
            if nargin > 0
                if isnumeric(com)
                    com = num2str(com);
                end
                com = upper(com);
                if ~strncmp(com,'COM',3)
                    com = ['COM' com];
                end
                obj.Ins = serial(com, 'BaudRate', 9600);
            else
                error('Please specify the COM port for the CS580.');
            end
        end
        
        function data = IDN(obj) %Indentify
        command = '*IDN?' ;
        fprintf(obj.Ins, command);
        data = fscanf(obj.Ins,'%c');
%         disp(data);
        end
        
        %% set functions: 
        
        function setCurrent(obj, curr) % for better understading of the gain options -> see the table below.
            % Calculate the order of magnitude of the input number
            orderOfMagnitude = log10(abs(curr));
            % Perform actions based on the order of magnitude
            if orderOfMagnitude > -1 || orderOfMagnitude < -13
                error('Not in the valid range')
            end
            if orderOfMagnitude <= -9
                order  = 0 ; 
            else
                if orderOfMagnitude >= -2
                    order = 8;
                else
                    order = floor(orderOfMagnitude + 10);
                end
            end
            command = ['GAIN ' num2str(order) ] ; % There is no problem to change the gain during higher or lower current, it's takes it down
            fprintf(obj.Ins, command);
            command = ['CURR ' num2str(curr)] ;
            fprintf(obj.Ins, command);
        end 
        
        function OutputON(obj)
            command = 'SOUT ON'; 
            fprintf(obj.Ins, command);
        end
        
        function OutputOFF(obj)
            command = 'SOUT OFF'; 
            fprintf(obj.Ins, command);
        end
        
        function setCurr_withoutFixingGain(obj, curr) % not recomended! we need to change the gain as well. use setCurr function instead.
        command = ['CURR ' num2str(curr)] ;
        fprintf(obj.Ins, command);
        end
        
        function setGain_withoutFixingCurrentAndVoltage(obj, gain) % not recomended! we need to change the current as well. use setCurr function instead.
        command = ['GAIN ' num2str(gain)] ; %notice that the gain need to be chosen from the next table:
        fprintf(obj.Ins, command);
        end
                % Values for gain option:
        % command name (z)   Gain            command name (z)     Gain
        % G1nA \(or) 0       1 nA/V          G100uA \ 5           100 uA/V
        % G10nA \  1         10 nA/V         G1mA \ 6             1 mA/V
        % G100nA \ 2         100 nA/V        G10mA \ 7            10 mA/V
        % G1uA  \ 3          1 uA/V          G50mA \ 8            50 mA/V
        % G10uA \ 4          10 uA/V
        
        function setVoltage_withoutFixingGainandCurrent(obj, volt) % not reconended! we need to change the gain as well. use setCurr function instead.
        % The default value is VOLT 10.0.
        % The allowable range for VOLT is always 0 to 50.
        command = ['VOLT ' num2str(volt)] ;
        fprintf(obj.Ins, command);
        end
        
        function InputON(obj)
            command = 'INPT ON'; 
            fprintf(obj.Ins, command); 
        end
        
        function InputOFF(obj)
            command = 'INPT OFF'; 
            fprintf(obj.Ins, command);
        end
        
        
        
        %% Get functions: 
                
        function response = getCurrent(obj)
        command = 'CURR?' ;
        fprintf(obj.Ins, command);
        response = fscanf(obj.Ins);
        response = str2num(response);
        end
        
        function response = getGain(obj)
        command = 'GAIN?' ;
        fprintf(obj.Ins, command);
        response = fscanf(obj.Ins);
        end
        
        function response = getVoltage(obj)
        command = 'VOLT?' ;
        fprintf(obj.Ins, command);
        response = fscanf(obj.Ins);
        response = str2num(response);
        end
        

        
    

    
    
    
    end
end
