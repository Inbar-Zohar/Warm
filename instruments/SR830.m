classdef SR830 < handle
    properties
        Ins
    end

    methods
        function obj = SR830(board,address)
            obj.Ins = gpib('NI', board, address, 'Timeout', 1.0);
        end

        function Clear(obj); fprintf(obj.Ins,'*CLS'); end % Clear
        function Reset(obj); fprintf(obj.Ins,'*RST'); end % Reset Configuration

        function response = robustQuery(obj, command, max_attempts, pause_time)
            if nargin < 3; max_attempts = 10; end
            if nargin < 4; pause_time = 0.01; end
            for attempt = 1:max_attempts
                flushinput(obj.Ins); % Clear input buffer
                flushoutput(obj.Ins); % Clear output buffer
                obj.Clear; % Clear device status
                response = query(obj.Ins, command);
            end
        end

        function data = IDN(obj); data = obj.robustQuery('*IDN?'); end % IDN
        function data = Device(obj); IDN = obj.IDN; data = IDN(27:31); end % Device model        

        %% Get current measurement results
        function data = ShowChannelScreen(obj,ch) %Show channel screen
            data = obj.robustQuery(sprintf('OUTR? %d',ch));
        end
        function data = ShowX(obj) %Show X
            data = str2double(obj.robustQuery('OUTP? 1'));
        end
        function data = ShowY(obj) %Show Y
            data = str2double(obj.robustQuery('OUTP? 2'));
        end
        function data = ShowR(obj) %Show R
            data = str2double(obj.robustQuery('OUTP? 3'));
        end
        function data = ShowTheta(obj) %Show omega
            data = str2double(obj.robustQuery('OUTP? 4'));
        end
        function data = ShowFreq(obj) %Show frequency
            data = str2double(obj.robustQuery('FREQ?'));
        end


        %% Set/Get Configuration       
        function setSensitivity(obj,num)
            if isnumeric(num)
                val = 26+round(log10(num)*3);
                if val < 0 || val > 27
                    error('Bad input value %g. Should be 2e-9<=x<=1',num)
                end
            else
                error('Non numeric input value');
            end
            fprintf(obj.Ins,sprintf('SENS %d',val)); 
        end
        function data = getSensitivity(obj)
            x = str2double(obj.robustQuery('SENS?'));
            data = round(10^((mod(x,3)+1)/3))*10^(floor(x/3)-9);
        end       

        function setTConst(obj,num)
            if isnumeric(num)
                val = round(log10(num)*2)+10;
                if val < 0 || val > 19
                    error('Bad input value %g. Should be 10e-6<=x<=30e3',num)
                end
            else
                error('Non numeric input value');
            end
            fprintf(obj.Ins,sprintf('OFLT %d',val)); 
        end
        function data = getTConst(obj)
            x = str2double(obj.robustQuery('OFLT?'));
            data = (1 + 2*mod(x,2)) .* 10.^(floor((x)/2) - 5);
        end       


        % NOTE: these commands can be used for current (i.e. not voltage) input, currently not supported
        function setSignalInputMode(obj,typ)
            if ischar(typ)
                if strcmpi(typ,'A')
                    val = 0;
                elseif strcmpi(typ,'A-B')
                    val = 1;
                else
                    error('String input should be "A" or "A-B", not %s',typ)
                end
            else
                val = double(logical(typ));
            end
            fprintf(obj.Ins,sprintf('ISRC %d',val));
        end
        function setSignalInputModeA(obj); obj.setSignalInputMode(0); end
        function setSignalInputModeAminusB(obj); obj.setSignalInputMode(1); end
        function typ = getSignalInputMode(obj)
            x = str2double(obj.robustQuery('ISRC?'));
            switch x; case 0; typ = 'A'; case 1; typ = 'A-B'; end
        end

        function setSignalInputCoupling(obj,typ)
            if ischar(typ)
                if strcmpi(typ,'AC')
                    val = 0;
                elseif strcmpi(typ,'DC')
                    val = 1;
                else
                    error('String input should be "AC" or "DC", not %s',typ)
                end
            else
                val = double(logical(typ));
            end
            fprintf(obj.Ins,sprintf('ICPL %d',val));
        end
        function setSignalInputCouplingAC(obj); obj.setSignalInputCoupling(0); end
        function setSignalInputCouplingDC(obj); obj.setSignalInputCoupling(1); end
        function typ = getSignalInputCoupling(obj)
            x = str2double(obj.robustQuery('ICPL?'));
            if x==1; typ = 'DC'; else; typ = 'AC'; end
        end

        function setSignalInputShield(obj,typ)
            if ischar(typ)
                if strcmpi(typ,'Float')
                    val = 0;
                elseif strcmpi(typ,'Ground')
                    val = 1;
                else
                    error('String input should be "Float" or "Ground", not %s',typ)
                end
            else
                val = double(logical(typ));
            end
            fprintf(obj.Ins,sprintf('IGND %d',val));
        end
        function setSignalInputShieldFloat(obj); obj.setSignalInputShield(0); end
        function setSignalInputShieldGround(obj); obj.setSignalInputShield(1); end
        function typ = getSignalInputShield(obj)
            x = str2double(obj.robustQuery('IGND?'));
            if ~x; typ = 'Float'; else; typ = 'Ground'; end
        end

        function setLineNotchFilter(obj,typ)
            if ischar(typ)
                if strcmpi(typ,'Off')
                    val = 0;
                elseif strcmpi(typ,'Line')
                    val = 1;
                elseif strcmpi(typ,'2xLine')
                    val = 2;
                elseif strcmpi(typ,'Both')
                    val = 3;
                else
                    error('String input should be "Off", "Line", "2xLine" or "Both", not %s',typ)
                end
            else
                val = typ;
            end
            fprintf(obj.Ins,sprintf('ILIN %d',val));
        end
        function setLineNotchFilterOff(obj); obj.setLineNotchFilter(0); end
        function setLineNotchFilterLine(obj); obj.setLineNotchFilter(1); end
        function setLineNotchFilter2xLine(obj); obj.setLineNotchFilter(2); end
        function setLineNotchFilterBoth(obj); obj.setLineNotchFilter(3); end

        function typ = getLineNotchFilter(obj)
            x = str2double(obj.robustQuery('ILIN?'));
            switch x
                case 0; typ = 'Off'; 
                case 1; typ = 'Line'; 
                case 2; typ = '2xLine'; 
                case 3; typ = 'Both'; 
            end
        end


        function setFilterOrder(obj,num)
            if isnumeric(num)
                val = round(num-1);
                if val < 0 || val > 3
                    error('Bad input value %g. Should be 1<=x<=4',num)
                end
            else
                error('Non numeric input value');
            end
            fprintf(obj.Ins,sprintf('OFSL %d',val)); 
        end
        function data = getFilterOrder(obj)
            x = str2double(obj.robustQuery('OFSL?'));
            data = x + 1;
        end       

        function setSyncFilter(obj,stt)
            val = double(logical(stt));
            fprintf(obj.Ins,sprintf('SYNC %d',val));
        end
        function setSyncFilterON(obj); obj.setSyncFilter(1); end
        function setSyncFilterOFF(obj); obj.setSyncFilter(0); end
        function stt = getSyncFilter(obj)
            stt = str2double(obj.robustQuery('SYNC?'));
        end

        % NOTE: these commands can be used for setting other than X,Y,R,Theta, currently not supported
        function setOutputCh(obj,ch,typ)
            if ch == 1
                if ischar(typ)
                    if strncmpi(typ,'X',1)
                        val = 0;
                    elseif strncmpi(typ,'R',1)
                        val = 1;
                    else
                        error('CH1 can be set to either X or R, requested %s', typ)
                    end
                else
                    val = double(logical(typ));
                end
            elseif ch == 2
                if ischar(typ)
                    if strcmpi(typ,'Y') || strcmpi(typ,'XY')
                        val = 0;
                    elseif strncmpi(typ,'TH',2) || strncmpi(typ,'RTH',3) 
                        val = 1;
                    else
                        error('CH2 can be set to either Y or THeta, requested %s', typ)
                    end
                else
                    val = double(logical(typ));
                end
            else
                error('ch must be 1 or 2, requested %g', ch)
            end
            fprintf(obj.Ins,sprintf('DDEF %d, %d, 0',ch, val));
        end               
        function setOutputCh1(obj,typ); obj.setOutputCh(1,typ); end
        function setOutputCh2(obj,typ); obj.setOutputCh(2,typ); end
        function setOutputChX(obj); obj.setOutputCh(1,'X'); end
        function setOutputChY(obj); obj.setOutputCh(2,'Y'); end
        function setOutputChR(obj); obj.setOutputCh(1,'R'); end
        function setOutputChTheta(obj); obj.setOutputCh(2,'TH'); end

        function data = getOutputCh(obj,ch)
            out = strsplit(obj.robustQuery(sprintf('DDEF? %d',ch)),',');
            val = str2double(out{1});
            if ch == 1
                if val == 0; data = 'X'; else; data = 'R'; end
            elseif ch == 2
                if val == 0; data = 'Y'; else; data = 'TH'; end
            end
        end
        function data = getOutputCh1(obj); data = getOutputCh(obj,1); end
        function data = getOutputCh2(obj); data = getOutputCh(obj,2); end

        % Set Output Expand
        function setOutputExpand(obj, ch, expand)
            if ischar(ch)
                if strcmpi(ch, 'X'); ch_ind = 1;  % SR830 uses 1=X, 2=Y, 3=R
                elseif strcmpi(ch, 'Y'); ch_ind = 2;
                elseif strcmpi(ch, 'R'); ch_ind = 3;
                else; error('ch must be X, Y or R, requested %s', ch);
                end
            elseif isnumeric(ch)
                if ismember(ch, [1 2 3]); ch_ind = ch;
                else; error('ch must be 1, 2 or 3, requested %g', ch);
                end
            else
                disp(ch)
                error('What the hell did you give as channel? disp(ch) is printed above')
            end
            if ischar(expand)
                if strcmpi(expand, 'OFF'); expand_ind = 0;
                elseif strcmpi(expand, '10') || strcmpi(expand, 'x10'); expand_ind = 1;
                elseif strcmpi(expand, '100') || strcmpi(expand, 'x100'); expand_ind = 2;
                else; error('expand must be OFF, x10 or x100, requested %s', expand);
                end
            elseif isnumeric(expand)
                if ismember(expand, [0 1 2]); expand_ind = expand;
                else; error('expand must be 0, 1 or 2, requested %g', expand);
                end
            else
                disp(expand)
                error('What the hell did you give as the output expand factor? disp(expand) is printed above')
            end
            % For SR830, we must set offset too; query current offset and preserve it
            current_offset = obj.getOutputOffsetValue(ch);  % Fetch current offset
            fprintf(obj.Ins, sprintf('OEXP %d, %.2f, %d', ch_ind, current_offset, expand_ind));
        end
        function setOutputExpandX(obj, expand); obj.setOutputExpand('X', expand); end
        function setOutputExpandY(obj, expand); obj.setOutputExpand('Y', expand); end
        function setOutputExpandR(obj, expand); obj.setOutputExpand('R', expand); end
        
        % Get Output Expand
        function data = getOutputExpand(obj, ch)
            if ischar(ch)
                if strcmpi(ch, 'X'); ch_ind = 1;
                elseif strcmpi(ch, 'Y'); ch_ind = 2;
                elseif strcmpi(ch, 'R'); ch_ind = 3;
                else; error('ch must be X, Y or R, requested %s', ch);
                end
            elseif isnumeric(ch)
                if ismember(ch, [1 2 3]); ch_ind = ch;
                else; error('ch must be 1, 2 or 3, requested %g', ch);
                end
            else
                disp(ch)
                error('What the hell did you give as channel? disp(ch) is printed above')
            end
            response = obj.robustQuery( sprintf('OEXP? %d', ch_ind));  % Returns "offset,expand"
            vals = split(strip(response), ',');  % Split into [offset, expand]
            expand_val = str2double(vals{2});
            switch expand_val
                case 0; data = 'OFF';
                case 1; data = 'x10';
                case 2; data = 'x100';
            end
        end
        function data = getOutputExpandX(obj); data = getOutputExpand(obj, 'X'); end
        function data = getOutputExpandY(obj); data = getOutputExpand(obj, 'Y'); end
        function data = getOutputExpandR(obj); data = getOutputExpand(obj, 'R'); end
        
        % Set Output Offset Value
        function setOutputOffsetValue(obj, ch, offset)
            if ischar(ch)
                if strcmpi(ch, 'X'); ch_ind = 1;
                elseif strcmpi(ch, 'Y'); ch_ind = 2;
                elseif strcmpi(ch, 'R'); ch_ind = 3;
                else; error('ch must be X, Y or R, requested %s', ch);
                end
            elseif isnumeric(ch)
                if ismember(ch, [1 2 3]); ch_ind = ch;
                else; error('ch must be 1, 2 or 3, requested %g', ch);
                end
            else
                disp(ch)
                error('What the hell did you give as channel? disp(ch) is printed above')
            end
            if abs(offset) > 105.00; error('offset must be between -105.00 and 105.00, requested %f', offset); end
            % For SR830, we must set expand too; query current expand and preserve it
            current_expand = obj.getOutputExpand(ch);
            if strcmpi(current_expand, 'OFF'); expand_ind = 0;
            elseif strcmpi(current_expand, 'x10'); expand_ind = 1;
            elseif strcmpi(current_expand, 'x100'); expand_ind = 2;
            end
            fprintf(obj.Ins, sprintf('OEXP %d, %.2f, %d', ch_ind, offset, expand_ind));
        end
        function setOutputOffsetValueX(obj, offset); obj.setOutputOffsetValue('X', offset); end
        function setOutputOffsetValueY(obj, offset); obj.setOutputOffsetValue('Y', offset); end
        function setOutputOffsetValueR(obj, offset); obj.setOutputOffsetValue('R', offset); end
        
        % Get Output Offset Value
        function data = getOutputOffsetValue(obj, ch)
            if ischar(ch)
                if strcmpi(ch, 'X'); ch_ind = 1;
                elseif strcmpi(ch, 'Y'); ch_ind = 2;
                elseif strcmpi(ch, 'R'); ch_ind = 3;
                else; error('ch must be X, Y or R, requested %s', ch);
                end
            elseif isnumeric(ch)
                if ismember(ch, [1 2 3]); ch_ind = ch;
                else; error('ch must be 1, 2 or 3, requested %g', ch);
                end
            else
                disp(ch)
                error('What the hell did you give as channel? disp(ch) is printed above')
            end
            response = obj.robustQuery( sprintf('OEXP? %d', ch_ind));  % Returns "offset,expand"
            vals = split(strip(response), ',');  % Split into [offset, expand]
            data = str2double(vals{1});  % Offset is first value
        end
        function data = getOutputOffsetValueX(obj); data = getOutputOffsetValue(obj, 'X'); end
        function data = getOutputOffsetValueY(obj); data = getOutputOffsetValue(obj, 'Y'); end
        function data = getOutputOffsetValueR(obj); data = getOutputOffsetValue(obj, 'R'); end


        function setPhase(obj,num)
            fprintf(obj.Ins,sprintf('PHAS %d',num));
        end
        function data = getPhase(obj)
            data = str2double(obj.robustQuery('PHAS?'));
        end

        function setRefSource(obj,typ)
            if ischar(typ)
                if strncmpi(typ,'internal',3)
                    val = 1;
                elseif strncmpi(typ,'external',3)
                    val = 0;
                else                    
                    error('String input should be "INT"/"EXT", not %s',typ)
                end
            else
                if ismember(typ,0:1)
                    val = typ;
                else
                    error('Numerical input should be 0 to 1, not %g',typ)
                end
            end
            fprintf(obj.Ins,sprintf('FMOD %d',val));
        end
        function data = getRefSource(obj)
            val = str2double(obj.robustQuery('FMOD?'));
            switch val
                case 0
                    data = 'EXT';
                case 1
                    data = 'INT';
            end
        end

        function setExtRefMode(obj,typ)
            if ischar(typ)
                if strncmpi(typ,'sine',3)
                    val = 0;
                elseif strncmpi(typ,'posTTL',3)
                    val = 1;
                elseif strncmpi(typ,'negTTL',3)
                    val = 2;
                else                    
                    error('String input should be "SIN"/"POS"/"NEG", not %s',typ)
                end
            else
                if ismember(typ,0:2)
                    val = typ;
                else
                    error('Numerical input should be 0 to 2, not %g',typ)
                end
            end
            fprintf(obj.Ins,sprintf('RSLP %d',val));
        end
        function data = getExtRefMode(obj)
            val = str2double(obj.robustQuery('RSLP?'));
            switch val
                case 0
                    data = 'SIN';
                case 1
                    data = 'POS';
                case 2
                    data = 'NEG';
            end
        end


        %% Misc
        function AutoPhase(obj); fprintf(obj.Ins,'APHS'); end


        %% Full configuration commands
        function [names,getters] = allConfigurationVariables(obj)
            full_array = {...
                'phase','Phase'
                'sensitivity','Sensitivity'
                'time_const','TConst'
                'signal_input_mode','SignalInputMode'
                'signal_input_coupling','SignalInputCoupling'
                'signal_input_shield','SignalInputShield'
                'filter_order','FilterOrder'
                'sync_filter','SyncFilter'
                'output_ch1','OutputCh1'
                'output_expand_x','OutputExpandX'
                'output_expand_r','OutputExpandR'   
                'output_offset_x_value','OutputOffsetValueX'
                'output_offset_r_value','OutputOffsetValueR'
                'output_ch2','OutputCh2'
                'output_expand_y','OutputExpandY'  
                'output_offset_y_value','OutputOffsetValueY'                
                'ref_source','RefSource'
                'ext_ref_mode','ExtRefMode'
                };

            names = full_array(:,1);
            getters = full_array(:,2);
        end

        function getter = findGetterByName(obj,name,typ)
            if ~exist('typ','var')
                typ = 'get';
            end
            typ = lower(typ(1:3));
            if ~(strcmp(typ,'get') || strcmp(typ,'set'))
                error('typ should define getter or setter, not %s', typ);
            end
                
            [names,getters] = obj.allConfigurationVariables;
            getter = [typ getters{strcmp(name,names)}];
        end

        function config = readConfiguration(obj) %Read configuration
            [names,getters] = obj.allConfigurationVariables;
            for ind=1:length(names)
                for bla=1:3 
                    % Something stupid this way comes
                    try; config.(names{ind}) = obj.(obj.findGetterByName(names{ind},'get')); end
                end
            end
        end

        function applyConfiguration(obj,config) %Read configuration
            fields = fieldnames(config);
            names = obj.allConfigurationVariables;
            for ind=1:length(fields)
                if ismember(fields{ind},names)
                    obj.(obj.findGetterByName(fields{ind},'set'))(config.(fields{ind}));
                end
            end
        end
    end
end