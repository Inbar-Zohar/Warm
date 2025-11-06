classdef TDS2014C
    properties(Transient)
        Ins
    end
    methods
        function obj = TDS2014C(address)
            obj.Ins = instrfind('Type', 'visa-usb', 'RsrcName', 'USB0::0x0699::0x03A4::C042901::0::INSTR', 'Tag', '');
            if isempty(obj.Ins)
                obj.Ins = visa('AGILENT',['USB0::' address '::0::INSTR']);
            else
                fclose(obj.Ins);
                obj.Ins = obj.Ins(1);
            end
            set(obj.Ins,'Timeout',2);
            set(obj.Ins,'InputBufferSize',2^18);
            set(obj.Ins,'OutputBufferSize',2^18);
        end
        
        function data = IDN(obj) %IDN
            fprintf(obj.Ins,'*IDN?');
            data = fscanf(obj.Ins);
        end
        
        function Reset(obj)
            fprintf(obj.Ins, '*RST');
        end
        
        function state = State(obj)
            fprintf(obj.Ins,'*OPC?');
            state = num2str(fscanf(obj.Ins));
        end
        
        function RunStop(obj,sw)
            if isnumeric(sw) || islogical(sw);
                if sw; sw = 'RUN'; else sw = 'STOP'; end;
            end
            fprintf(obj.Ins, ['ACQ:STATE ',sw]);
        end
        
        function Run(obj)
            fprintf(obj.Ins, 'ACQ:STATE 1');
        end
        
        function Single(obj)
            fprintf(obj.Ins,'ACQuire:STOPAfter SEQ');
            fprintf(obj.Ins,'ACQ:STATE ON');
        end
        
        function TrigSource(obj,trig) %trig='CH1','CH2','CH3','CH4','EXT','LINE'
            fprintf(obj.Ins, ['TRIG:MAI:EDGE:SOU ',trig]);
        end
        
        function TrigMode(obj,trigM) %trigM= 'AUTO' | 'NORMal'
            fprintf(obj.Ins, ['TRIG:MAI:MODE ',trigM]);
        end
        
        function ForceTrig(obj); fprintf(obj.Ins,'TRIGger FORCe'); end
        
        function S = GetSensitivity(obj,ch) %Ch=[1,3,4]
            S = zeros(size(ch));
            for Ind = 1:length(ch)
                fprintf(obj.Ins,['CH' num2str(ch(Ind)) ':VOLTS?']);
                S(Ind) = str2double(fscanf(obj.Ins));
            end
        end
        
        function SetSensitivity(obj,ch,snstvty) %Ch=[1,3,4] snstvty=[0.1,2.5,0.6]
            for Ind = 1:length(ch)
                fprintf(obj.Ins,...
                    ['CH' num2str(ch(Ind)) ':VOLTS ' num2str(snstvty(Ind))]);
            end
        end
        
        function setTscale(obj,dt)
            fprintf(obj.Ins,['HOR:MAIN:SCALE ' num2str(dt,'%E')]);
        end
        
        function setVscale(obj,ch,Vscale)
            fprintf(obj.Ins, sprintf([':CH',num2str(ch),':SCAL ',num2str(Vscale)]));
        end

        function scale = getTscale(obj)
            scale = str2double(query(obj.Ins,'HOR:MAIN:SCALE?'));
        end
        
        function scale = getVscale(obj,ch)
            scale = nan(size(ch));
            for Ind = 1:length(ch)
            scale(Ind) = str2double(query(obj.Ins,[':CH',num2str(ch(Ind)),':SCAL?']));
            end
        end
        
        function [t,v] = Read(obj,ch) % ch=[1,3]
            fprintf(obj.Ins,['DATa:SOUrce CH' num2str(ch(1))]); obj.OPC;
            S = query(obj.Ins,'WFMP?'); obj.OPC;
            Ind = regexp(S,'points');
            N_points = str2double(S((find(S(1:Ind)==',',1,'last')+1):(Ind-1)));
%             RunStop(obj,'STOP'); %RUN & STOP
            
            % time vector:
            dt = str2double(query(obj.Ins,'WFMPre:XINcr?')); % obj.OPC;
            t0 = str2double(query(obj.Ins,'WFMPre:XZEro?')); % obj.OPC;
            t = t0+dt*((1:N_points)-1);
            
            % channels vectors
            if length(ch)==1; v = nan(size(t)); else v = nan(length(ch),N_points); end
            for Ind = 1:length(ch)
                if ~strcmpi(query(obj.Ins,'DATa:SOUrce?'),sprintf('CH%d\n',ch(Ind)))
                    fprintf(obj.Ins,['DATa:SOUrce CH' num2str(ch(Ind))]);
                end
                if ~strcmpi(query(obj.Ins,'DATa:ENC?'),sprintf('ASCII\n'))
                    fprintf(obj.Ins,'DATa:ENCdg ASCIi'); pause(2);
                end
                pause(1);
                fprintf(obj.Ins,'CURVe?');
                pause(4);
                data = fgetl(obj.Ins);
                if data(end)==10; data(end) = []; end
                data = char(data);
                data = str2double(strsplit(data,','));
                v0 = str2double(query(obj.Ins,'WFMPre:Yzero?'));
                vScale = str2double(query(obj.Ins,'WFMPre:YMUlt?'));
                fprintf(obj.Ins,'WFMPre:YOFf?');
                vOffset = str2double(fscanf(obj.Ins));
%                 if length(ch)==1; ch = 1; end;
                v(Ind,:) = v0 + vScale*(data-vOffset);
            end
            
            %           RunStop(obj,'RUN'); %RUN & STOP
        end
        
        function state = OPC(obj)
            fprintf(obj.Ins,'*OPC?'); state = fscanf(obj.Ins);
        end
        
        function err = readError(obj,dispFlag)
            err = query(obj.Ins,'ERRLOG:NEXT?');
            if exist('dispFlag','var') && dispFlag; disp(err); end
        end
        
        function offset = getVoffset(obj,ch)
            offsetPos = nan(size(ch));
            for Ind = 1:length(ch)
            offsetPos(Ind) = str2double(query(obj.Ins,[':CH',num2str(ch(Ind)),':POS?']));
            end
            scale = getVscale(obj,ch);
            offset = -offsetPos.*scale;
        end
        
        function setVoffset(obj,ch,offset)
            scale = getVscale(obj,ch);
            offsetPos = -offset./scale;
            for Ind = 1:length(ch)
            fprintf(obj.Ins,[':CH',num2str(ch(Ind)),':POS ' num2str(offsetPos(Ind))]);
            end
        end
        
        function R = getVrange(obj,ch)
            if ~exist('ch','var'); ch = 1:4; end;
            V0 = getVoffset(obj,ch); dV = getVscale(obj,ch);
            R = [V0(:)-5*dV(:),V0(:)+5*dV(:)];
%             if length(ch)==1; R = R(ch,:); end;
        end
        
        function setVrange(obj,ch,range)
            if ~exist('ch','var') || ~isnumeric(ch) || isnan(ch); ch = 1:4; end;
            for Ind = 1:length(ch)
                obj.setVscale(ch(Ind),diff(minmax(range(Ind,:)))/10);
                obj.setVoffset(ch(Ind),mean(range(Ind,:)));
            end
        end
        
        function setVrange4signal(obj,ch,y,scale_factor)
            if ~exist('scale_factor','var'); scale_factor = 1/1.3; end;
            range = minmax(y(:).'); range = 1/scale_factor*(range-mean(range))+mean(range);
            if diff(range)<0.013; range = mean(range)+0.010*[-1,1];
            elseif diff(range)>0.013 && diff(range)<0.020; range = mean(range)+0.020*[-1,1];
            end;
            if ~exist('ch','var') || ~isnumeric(ch) || isnan(ch); ch = 1:4; end;
            obj.setVrange(ch,range)
        end
        
    end
end