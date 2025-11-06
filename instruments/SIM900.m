classdef SIM900
    % TODO: 
    % add verifications that the card in the called port is the correct one

    properties
        Ins
    end

    methods
        function obj = SIM900(address)
            if ischar(address)
                if contains(address, 'COM', 'IgnoreCase', true)
                    comport = address;
                else
                    comport = sprintf('COM%s', address);
                end
            elseif isnumeric(address)
                comport = sprintf('COM%d', address);
            end
            obj.Ins = serial(comport, 'BaudRate', 9600);

            obj.Ins.Timeout = 4;
        end

        function data = IDN(obj,port) %IDN
            if ~exist('port','var') || isempty(port) || isnan(port)
                fprintf(obj.Ins,'*IDN?');
                data = fscanf(obj.Ins);
            else
                fprintf(obj.Ins,['conn ' num2str(port) ' , ''!' num2str(port) 'XYZ''']);
                fprintf(obj.Ins,'*IDN?');
                data = fscanf(obj.Ins);
                fprintf(obj.Ins,['!' num2str(port) 'XYZ']);
            end
        end

        function printToPort(obj,port,cmd)
            fprintf(obj.Ins,['conn ' num2str(port) ' ,''!' num2str(port) 'XYZ''']);
            fprintf(obj.Ins,cmd);
            fprintf(obj.Ins,['!' num2str(port) 'XYZ']);
        end

        function data = queryPort(obj,port,cmd)
            fprintf(obj.Ins,['conn ' num2str(port) ' ,''!' num2str(port) 'XYZ''']);
            data = str2double(query(obj.Ins,cmd));
            fprintf(obj.Ins,['!' num2str(port) 'XYZ']);
        end


        %% get functions

        function data = getPropState(obj,port) %Proportional Gain
            data = obj.queryPort(port,'PCTL?');
        end

        function data = getPropGain(obj,port) %Proportional Gain
            data = obj.queryPort(port,'GAIN?');
        end

        function data = getIntState(obj,port) %Integral Gain
            data = obj.queryPort(port,'ICTL?');
        end

        function data = getIntGain(obj,port) %Integral Gain
            data = obj.queryPort(port,'INTG?');
        end

        function data = getDerState(obj,port) %Derivative Gain
            data = obj.queryPort(port,'DCTL?');
        end

        function data = getDerGain(obj,port) %Derivative Gain
            data = obj.queryPort(port,'DERV?');
        end

        function data = getOutputMode(obj,port) %Proportional Gain
            data = obj.queryPort(port,'AMAN?');
        end

        % function data = getOffset(obj,port) % Offset value
        %   fprintf(obj.Ins,['conn ' num2str(port) ' ,''!' num2str(port) 'XYZ''']);
        %   fprintf(obj.Ins,'OFST?');
        %   data=str2double(fscanf(obj.Ins));
        %   fprintf(obj.Ins,['!' num2str(port) 'XYZ']);
        % end
        %
        % function data = getUpLIM(obj,port) %Controller Upper limit output [Volts]
        %   fprintf(obj.Ins,['conn ' num2str(port) ' ,''!' num2str(port) 'XYZ''']);
        %   fprintf(obj.Ins,'ULIM?');
        %   data=str2double(fscanf(obj.Ins));
        %   fprintf(obj.Ins,['!' num2str(port) 'XYZ']);
        % end
        %
        %
        % function data = getLowLIM(obj,port) %Controller Lower limit output [Volts]
        %   fprintf(obj.Ins,['conn ' num2str(port) ' ,''!' num2str(port) 'XYZ''']);
        %   fprintf(obj.Ins,'LLIM?');
        %   data=str2double(fscanf(obj.Ins));
        %   fprintf(obj.Ins,['!' num2str(port) 'XYZ']);
        % end


        %% Set functions
        function setPIDState(obj,port,typ,stt) %Turn on (on) /off (off) Proportional Gain
            if ischar(stt)
                if strcmpi(stt,'OFF')
                    val = 0;
                elseif strcmpi(stt,'ON')
                    val = 1;
                else
                    error('String input should be "ON" or "OFF", not %s',stt)
                end
            else
                val = double(logical(stt));
            end
            if ~ismember(typ,{'P','I','D'})
                error('type should be "P","I" or "D", not %s',typ)
            end
            obj.printToPort(port,sprintf('%sCTL %d',typ,val));
        end

        function setPropState(obj,port,stt); setPIDState(obj,port,'P',stt); end
        function setPropStateON(obj,port); obj.setPropState(port,'ON'); end
        function setPropStateOFF(obj,port); obj.setPropState(port,'OFF'); end

        function setPropGain(obj,port,gain,only_neg) %set Proportional Gain
            if ~exist('only_neg','var'); only_neg = 0; end
            if only_neg && gain >= 0; error("By default only negative gain allowed"); end
            obj.printToPort(port,['GAIN ' num2str(gain)]);
        end

        function setIntState(obj,port,stt); setPIDState(obj,port,'I',stt); end
        function setIntStateON(obj,port); obj.setIntState(port,'ON'); end
        function setIntStateOFF(obj,port); obj.setIntState(port,'OFF'); end

        function setIntGain(obj,port,gain) %set Integral Gain
            obj.printToPort(port,['INTG' num2str(gain)]);
        end

        function setDerState(obj,port,stt); setPIDState(obj,port,'D',stt); end
        function setDerStateON(obj,port); obj.setIntState(port,'ON'); end
        function setDerStateOFF(obj,port); obj.setIntState(port,'OFF'); end

        function setDerGain(obj,port,gain) %set Derivative Gain
            obj.printToPort(port,['DERV' num2str(gain)]);
        end

        function setOutputMode(obj,port,stt) % Choose between Manual/PID control
            if ischar(stt)
                if strcmpi(stt,'MAN')
                    val = 0;
                elseif strcmpi(stt,'PID')
                    val = 1;
                else
                    error('String input should be "MAN" or "PID", not %s',stt)
                end
            else
                val = double(logical(stt));
            end
            obj.printToPort(port,sprintf('AMAN %d',val));
        end        
        function setOutputModeMAN(obj,port); obj.setOutputMode(port,'MAN'); end
        function setOutputModePID(obj,port); obj.setOutputMode(port,'PID'); end        
            
        % function setOffset(obj,port,ofst) %set Offset value
        %     % The allowed range is  +-10.000 [V]
        %     fprintf(obj.Ins,['conn ' num2str(port) ' ,''!' num2str(port) 'XYZ''']);
        %     fprintf(obj.Ins,['OFST' num2str(ofst)]);
        %     fprintf(obj.Ins,'OFST?');
        %     data=str2double(fscanf(obj.Ins));
        %     fprintf(obj.Ins,['!' num2str(port) 'XYZ']);
        % end
        % 
        % function data = setInternalSetpoint(obj,port,setp) %set setpoint value
        %     % The allowed range is  +-10.000 [V]
        %     fprintf(obj.Ins,['conn ' num2str(port) ' ,''!' num2str(port) 'XYZ''']);
        %     fprintf(obj.Ins,['SETP' num2str(setp)]);
        %     fprintf(obj.Ins,'SETP?');
        %     data=str2double(fscanf(obj.Ins));
        %     fprintf(obj.Ins,['!' num2str(port) 'XYZ']);
        % end
        % 
        % function data = setRampingRate(obj,port,rate) %set ramping rate
        %     % The allowed range is  1E-3 <= rate <= 1E4
        %     fprintf(obj.Ins,['conn ' num2str(port) ' ,''!' num2str(port) 'XYZ''']);
        %     fprintf(obj.Ins,['RATE' num2str(rate)]);
        %     fprintf(obj.Ins,'RATE?');
        %     data=str2double(fscanf(obj.Ins));
        %     fprintf(obj.Ins,['!' num2str(port) 'XYZ']);
        % end
        % 
        % function data = setManualOutput(obj,port,outp) %set the manual output value
        %     % The allowed range is  +-10.000 [V]
        %     fprintf(obj.Ins,['conn ' num2str(port) ' ,''!' num2str(port) 'XYZ''']);
        %     fprintf(obj.Ins,['MOUT' num2str(outp)]);
        %     fprintf(obj.Ins,'MOUT?');
        %     data=str2double(fscanf(obj.Ins));
        %     fprintf(obj.Ins,['!' num2str(port) 'XYZ']);
        % end
        % 
        % function data = setUpLIM(obj,port,lim) %Controller Upper limit output [Volts]
        %     fprintf(obj.Ins,['conn ' num2str(port) ' ,''!' num2str(port) 'XYZ''']);
        %     fprintf(obj.Ins,['ULIM ' num2str(lim)]);
        %     fprintf(obj.Ins,'ULIM?');
        %     data=str2double(fscanf(obj.Ins));
        %     fprintf(obj.Ins,['!' num2str(port) 'XYZ']);
        % end
        % 
        % 
        % function data = setLowLIM(obj,port,lim) %Controller Lower limit output [Volts]
        %     fprintf(obj.Ins,['conn ' num2str(port) ' ,''!' num2str(port) 'XYZ''']);
        %     fprintf(obj.Ins,['LLIM ' num2str(lim)]);
        %     fprintf(obj.Ins,'LLIM?');
        %     data=str2double(fscanf(obj.Ins));
        %     fprintf(obj.Ins,['!' num2str(port) 'XYZ']);
        % end
        % 
        % 
        % function data = setSetPointMode(obj,port,state) % Choose between Internal/External Setpoint
        %     % state should be number or single letter -> INT(i)(0) / EXT (e)(1)
        %     fprintf(obj.Ins,['conn ' num2str(port) ' ,''!' num2str(port) 'XYZ''']);
        %     if exist('state','var') && ~isempty(state) && ~isnan(state)
        %         if ischar(state)
        %             switch lower(state); case 'i'; state = 0; case 'e'; state = 1; otherwise; error('ERROR! Setpoint status is either 0 or 1'); end
        %         end
        %         state = double(state);
        %         fprintf(obj.Ins,['INPT ' num2str(state)]);
        %     end
        %     data = str2double(query(obj.Ins,'AMAN?')); %
        %     fprintf(obj.Ins,['!' num2str(port) 'XYZ']); %is not relaible
        % end

        %% Full configuration commands
        function [names,getters] = allConfigurationVariables(obj)
            full_array = {...
                'prop_state','PropState'
                'prop_gain','PropGain'
                'int_state','IntState'
                'int_gain','IntGain'
                'der_state','DerState'
                'der_gain','DerGain'     
                'output_mode','OutputMode'
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

        function config = readConfiguration(obj,port) %Read configuration
            [names,getters] = obj.allConfigurationVariables;            
            for ind=1:length(names)
                config.(names{ind}) = obj.(obj.findGetterByName(names{ind},'get'))(port);
            end
        end

        function applyConfiguration(obj,port,config) %Read configuration
            fields = fieldnames(config);
            names = obj.allConfigurationVariables;
            for ind=1:length(fields)
                if ismember(fields{ind},names)
                    obj.(obj.findGetterByName(fields{ind},'set'))(port,config.(fields{ind}));
                end
            end
        end




        %% Controller configuration


        %       function data = getPIDA_IntGain(obj);  data = obj.ShowIntGain(obj.PIDAport);  end % returns Xe-129 Integral Gain
        %       function data = getPIDB_IntGain(obj);  data = obj.ShowIntGain(obj.PIDBport);  end % returns Xe-131 Integral Gain
        %       function data = getPIDA_PropGain(obj); data = obj.ShowPropGain(obj.PIDAport); end % returns Xe-129 Proportional Gain
        %       function data = getPIDB_PropGain(obj); data = obj.ShowPropGain(obj.PIDBport); end % returns Xe-131 Proportional Gain
        %       function data = getPIDA_DerGain(obj);  data = obj.ShowDerGain(obj.PIDAport);  end % returns Xe-129 Derivative Gain
        %       function data = getPIDB_DerGain(obj);  data = obj.ShowDerGain(obj.PIDBport);  end % returns Xe-131 Derivative Gain
        %
        %       function data = setPIDA_IntGain(obj,g);  data = obj.SetIntGain(obj.PIDAport,g);  end % returns Xe-129 Integral Gain
        %       function data = setPIDB_IntGain(obj,g);  data = obj.SetIntGain(obj.PIDBport,g);  end % returns Xe-131 Integral Gain
        %       function data = setPIDA_PropGain(obj,g); data = obj.SetPropGain(obj.PIDAport,g); end % returns Xe-129 Proportional Gain
        %       function data = setPIDB_PropGain(obj,g); data = obj.SetPropGain(obj.PIDBport,g); end % returns Xe-131 Proportional Gain
        %       function data = setPIDA_DerGain(obj,g);  data = obj.SetDerGain(obj.PIDAport,g);  end % returns Xe-129 Derivative Gain
        %       function data = setPIDB_DerGain(obj,g);  data = obj.SetDerGain(obj.PIDBport,g);  end % returns Xe-131 Derivative Gain


        %       function data = PIDAOutput(obj,state)
        %           if ~exist('state','var'); state = nan; end
        %           data = obj.Output(obj.PIDAport,state);
        %       end
        %       function data = PIDBOutput(obj,state)
        %           if ~exist('state','var'); state = nan; end
        %           data = obj.Output(obj.PIDBport,state);
        %       end

    end
end


% PCTL(?) z 3 ? 10 Proportional action ON/OFF
% ICTL(?) z 3 ? 10 Integral action ON/OFF
% DCTL(?) z 3 ? 10 Derivative action ON/OFF
% OCTL(?) z 3 ? 10 Offset ON/OFF
% GAIN(?) {f} 3 ? 10 Proportional Gain
% APOL(?) z 3 ? 11 Controller Polarity
% INTG(?) {f} 3 ? 11 Integral Gain
% DERV(?) {f} 3 ? 11 Derivative Gain
% OFST(?) {f} 3 ? 11 Output Offset
% Controller Configuration
% AMAN(?) z 3 ? 12 Output (Manual Output/PID Control)
% INPT(?) z 3 ? 12 Input (Internal/External Setpoint)
% SETP(?) {f} 3 ? 12 New setpoint
% RAMP(?) z 3 ? 12 Internal setpoint ramping ON/OFF
% RATE(?) {f} 3 ? 12 Setpoint ramping Rate
% RMPS? 3 ? 13 Setpoint ramping status
% STRT z 3 ? 13 Pause or continue ramping
% MOUT(?) {f} 3 ? 13 Manual Output
% ULIM(?) {f} 3 ? 13 Upper Output Limit
% LLIM(?) {f} 3 ? 14 Lower Output Limit