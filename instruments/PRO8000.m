% Nitzan Haviv 09.08.23

%Example code: in the same directory -> TEST_pro8000.m
% Configuration:
%  pro8000_1 = PRO8000('COMx', Active_slots_array ); 
%  fopen(pro8000_1.Ins);

classdef PRO8000
    
    % Make sure to change the configuration manually to RS-232 from the
    % default IEEE488
    
    properties
        Ins % The COM port used to connect the device
        Active_slots
    end
    
    methods
        % Constructor
        function obj = PRO8000(address, Active_slots)
            switch nargin
                case 0
                    error('Please specify the COM port for the PRO8000.');
                case 1
                    if ischar(address)
                        if contains(address, 'COM', 'IgnoreCase', true)
                            obj.Ins = serial( address, 'BaudRate', 19200);
                        end
                    end
                    if isnumeric(address)
                        obj.Ins = serial(sprintf('COM%d', Ins),'BaudRate', 19200);
                    end
                    obj.Active_slots = 1:8; % Default value

                case 2
                    if ischar(address)
                        if contains(address, 'COM', 'IgnoreCase', true)
                               obj.Ins = serial( address, 'BaudRate', 19200);
                        end
                    end
                    if isnumeric(address)
                               obj.Ins = serial(sprintf('COM%d', address),'BaudRate', 19200);
                    end
                    obj.Active_slots = Active_slots;
            end
            if nargin >= 3 
                error('To many inputs arguments')
            end
        end
        
        function data = IDN(obj) %Indentify
        command = '*IDN?' ;
        fprintf(obj.Ins, command);
        data = fscanf(obj.Ins,'%c');
%       disp(data);
        end
        
        %% set functions: 
        
        function setCurrent(obj, slot, curr) %in [A] 
            command = [':SLOT ' num2str(slot)] ;
            fprintf(obj.Ins, command);
            % Assuring the input curr is below the defined I_max
            command =  ":LIMC:SET?";
            fprintf(obj.Ins, command);
            limit = fscanf(obj.Ins); % has the shape: ':LIMC:SET 1.28999993E-01'

            % For some reason sometimes the response to the limit query is
            % the current setpoint. This is a hack that solves it, probably
            % something better should be done.      
            max_tries = 10;
            this_try = 1;
            while this_try<max_tries && ~contains(limit,'SET')  
                fprintf(obj.Ins, command);
                limit = fscanf(obj.Ins); % has the shape: ':LIMC:SET 1.28999993E-01'                
            end
            % If the hack didn't work
            if ~contains(limit,'SET')
                error('Cannot get Imax from device for some reason')
            end
                                  
            % Separating I_max from the string
            pattern = '-?\d+(\.\d+)?(?:[Ee][+\-]?\d+)?';
            % Extract the number using regular expression
            matches = regexp(limit, pattern, 'match');
            if ~isempty(matches)
                % Convert the matched string to a double
                I_max = str2double(matches{1});
            else
                error('No valid number found in the response - error in code!');
            end
            if curr <= I_max
                        command = [':ILD:SET ' num2str(curr)] ;
                        fprintf(obj.Ins, command);
            else
                error('Assigned current is larger than the defined I_max')
            end
%             warning( 'Set point calibration has been made. Check the setpoint ')
        end 
        
        function OutputON(obj, slot) 
            command = [':SLOT ' num2str(slot)];
            fprintf(obj.Ins, command);
            % Assuring the input curr is below the defined I_max
            command =  ":LIMC:SET?";
            fprintf(obj.Ins, command);
            response = fscanf(obj.Ins); % has the shape: ':LIMC:SET 1.28999993E-01'
            % Separating I_max from the string
            pattern = '-?\d+(\.\d+)?(?:[Ee][+\-]?\d+)?';
            % Extract the number using regular expression
            matches = regexp(response, pattern, 'match');
            if ~isempty(matches)
                % Convert the matched string to a double
                I_max = str2double(matches{1});
            else
                error('No valid number found in the response - error in code!');
            end
            command =  ':ILD:SET?';
            fprintf(obj.Ins, command);
            response = fscanf(obj.Ins); % has the shape: ':LIMC:SET 1.28999993E-01'
            % Separating I_max from the string
            pattern = '-?\d+(\.\d+)?(?:[Ee][+\-]?\d+)?';
            % Extract the number using regular expression
            matches = regexp(response, pattern, 'match');
            if ~isempty(matches)
                % Convert the matched string to a double
                I_set = str2double(matches{1});
            else
                error('No valid number found in the response - error in code!');
            end
            if I_set <= I_max
                command = ':LASER ON';
                fprintf(obj.Ins, command);
            else
                error('Assigned current is larger than the defined I_max')
            end
        end
        
        function OutputON_All(obj) % Unconnected slots will stay off
            for slot = obj.Active_slots
                OutputON(obj,slot)
            end
        end
        
        function OutputOFF(obj, slot)
            command = [':SLOT ' num2str(slot)];
            fprintf(obj.Ins, command);
            command = ':LASER OFF';
            fprintf(obj.Ins, command);
        end
        
        function OutputOFF_All(obj)
            for slot = obj.Active_slots
              OutputOFF(obj, slot)  
            end
        end
        
%         % set the sensitivity calibration factor of the monitor diode [A/W]
%         function setCalibrationFactor(obj,slot, A_W)
%             command = [':SLOT ' num2str(slot)];
%             fprintf(obj.Ins, command);
%             command = [':CALPD:SET ' num2str(A_W)];
%             fprintf(obj.Ins, command);
%         end
        
        % set I_max for the specific slot
        function setI_max(obj,slot, I_max)
            command = [':SLOT ' num2str(slot)];
            fprintf(obj.Ins, command);
            % Assuring the input I_max is below I_Lim (hardware)
            command =  ":LIMCP:ACT?";
            fprintf(obj.Ins, command);
            response = fscanf(obj.Ins); % has the shape: ':LIMC:SET 1.28999993E-01'
            % Separating I_max from the string
            pattern = '-?\d+(\.\d+)?(?:[Ee][+\-]?\d+)?';
            % Extract the number using regular expression
            matches = regexp(response, pattern, 'match');
            if ~isempty(matches)
                % Convert the matched string to a double
                I_Lim = str2double(matches{1});
            else
                error('No valid number found in the response - error in code!');
            end
            if I_max < I_Lim
                command = [':LIMC:SET ' num2str(I_max)];
                fprintf(obj.Ins, command);
            else
                error('Assigned I_max is larger than I_Lim (hardware)')  
            end
        end
        
        % set constant current mode
        function setModeConstCurrent(obj ,slot)
            command = [':SLOT ' num2str(slot)];
            fprintf(obj.Ins, command);
            command = ":MODE CC";
            fprintf(obj.Ins, command);
        end
        
        % set constant power mode
        function setModeConstPower(obj ,slot)
            command = [':SLOT ' num2str(slot)];
            fprintf(obj.Ins, command);
            command = ":MODE CP";
            fprintf(obj.Ins, command);
        end
        
        % set the Laser diode polarity (anode / cathode grounded)
        function setAG_anode_grounded(obj ,slot)
%             command = [':SLOT ' num2str(slot)];
%             fprintf(obj.Ins, command);
%             command = ":LASER?";
%             fprintf(obj.Ins, command);
%             response = fscanf(obj.Ins);
            response = obj.getSlotStatus(slot);
            if contains(response, 'ON')
                error('Cannot alter diode polarity while the channel is ON');
            else
                command = ":LDPOL AG";
                fprintf(obj.Ins, command);
            end
        end
        
        function setCG_cathode_grounded(obj ,slot)
            response = obj.getSlotStatus(slot);
            if contains(response, 'ON')
                error('Cannot alter diode polarity while the channel is ON');
            else
                command = ":LDPOL CG";
                fprintf(obj.Ins, command);
            end
        end
        
        
        
        %% Get functions: 
        
        function response = getSlotStatus(obj,slot)
            command = [':SLOT ' num2str(slot)];
            fprintf(obj.Ins, command);
            command = ":LASER?";
            fprintf(obj.Ins, command);
            response = fscanf(obj.Ins);
%             disp(response);
        end
        
        
        % get the current set (ILD) in the current active slot
        function setPoint = getCurrentSetpoint(obj, slot)
        command = [':SLOT ' num2str(slot)];
        fprintf(obj.Ins, command);
        command = ':ILD:SET?';
        fprintf(obj.Ins, command);
        response = fscanf(obj.Ins);
        pattern = '-?\d+(\.\d+)?(?:[Ee][+\-]?\d+)?';
        matches = regexp(response, pattern, 'match');
        if ~isempty(matches)
            setPoint = str2double(matches{1});
        else
            error('No valid number found in the response - error in code!');
        end
        end
        
        % gets the actual laser diode current
        function current = getCurrent(obj, slot)
        command = [':SLOT ' num2str(slot)];
        fprintf(obj.Ins, command);
        command = ':ILD:ACT?';
        fprintf(obj.Ins, command);
        response = fscanf(obj.Ins);
        pattern = '-?\d+(\.\d+)?(?:[Ee][+\-]?\d+)?';
        matches = regexp(response, pattern, 'match');
        if ~isempty(matches)
            current = str2double(matches{1});
        else
            error('No valid number found in the response - error in code!');
        end
        end
        
%         % gets the sensitivity calibration factor of the monitor diode[A/W]
%         function response = getCalibrationFactor(obj,slot)
%         command = [':SLOT ' num2str(slot)];
%         fprintf(obj.Ins, command);
%         command = ':CALPD:SET?';
%         fprintf(obj.Ins, command);
%         response = fscanf(obj.Ins);
% %         disp(response);
%         end
        
        % gets I_max for the specific slot
        function I_max = getI_max(obj,slot)
        command = [':SLOT ' num2str(slot)];
        fprintf(obj.Ins, command);
        command =  ":LIMC:SET?";
        fprintf(obj.Ins, command);
        response = fscanf(obj.Ins);
        pattern = '-?\d+(\.\d+)?(?:[Ee][+\-]?\d+)?';
        matches = regexp(response, pattern, 'match');
        if ~isempty(matches)
            I_max = str2double(matches{1});
        else
            error('No valid number found in the response - error in code!');
        end
        end

        % gets the monitor diode current limit for ILIM - ADC = FFFF
        % Actual Hardware limit - for the specific slot
        % Can be physically modified with a screwdriver (under the serial
        % connector)
        function I_Lim = getHardwareI_Lim(obj,slot)
        command = [':SLOT ' num2str(slot)];
        fprintf(obj.Ins, command);
        command =  ":LIMCP:ACT?";
        fprintf(obj.Ins, command);
        response = fscanf(obj.Ins); 
        pattern = '-?\d+(\.\d+)?(?:[Ee][+\-]?\d+)?';
        matches = regexp(response, pattern, 'match');
        if ~isempty(matches)
            I_Lim = str2double(matches{1});
        else
            error('No valid number found in the response - error in code!');
        end
        end   
        
        
        
%         % gets the monitor diode current limit for ILIM - ADC = FFFF
%         % same for all slots
%         function response = getLimCurrent(obj,slot)
%         command = [':SLOT ' num2str(slot)];
%         fprintf(obj.Ins, command);
%         command =  ":LIMC:MAX_W?";
%         fprintf(obj.Ins, command);
%         response = fscanf(obj.Ins);
%         end
        

    end
end
