classdef DSOX3000
% Updates:
% 08/01/19 - YAHEL:
%	1- New option to read scope inital parameters from xlsx in the format of 'param2scope_yahel'
%	2- Adress 7 for my scope (with 'ni' vendor)
% 09/01/19 - ROY:
%   1- Edited READ function to be self-contained
%   2- Fixed YAHEL's scope #7 bug
%   3- Changed scope adresses table to SWITCH (from IF...ELSEIF...END)
%   4- Changed function "complete" to "Complete"
% 13/01/19 - YAHEL: Changed N_average>1 (instead of >0), line 124
% 07/05/19 - ROY: Added functions setVrange
% 28/08/19 - ROY: Added "high-resolution" to obj.Init, line 294
% 03/09/19 - ROY: Added getTdelay function, line 197
% 04/10/19 - ROY: Added functions setVrange4signal, line 134
% 19/10/19 - ROY: Added optional scale factor to setVrange4signal, line 134
% 02/12/19 - ROY: Added query option to TrigSource, line 229
% 03/12/19 - ROY: Added function Probe and added it to init_cahnnels and to Set, line 406
% 26/01/20 - ROY: Changed myBin to de2bi in Status(obj), line 107
% 26/01/20 - ROY: Added scopes 10,11, line 43
% 02/02/20 - ROY: Added function setChDisp, line 420
% 04/02/20 - ROY: Added function highres, line 423
% 04/02/20 - ROY: Changed AverageAQN for N_average=1 from NORMal to highres acquisition, line 163
% 26/03/20 - ROY: Added function getTdelay, line 232
% 31/03/20 - ROY: Added function getChImpedance, line 442
% 17/04/20 - ROY: Added function getChCoupling, line 450
% 29/06/20 - Avraham: Added case 14 to the scope constructor
% 08/07/20 - Gil: updated to RS functions from dor b code
% 28/05/21 - ROY: added default timeout to "readyToRead", line 118
% 02/09/24 - Kostya: added default setup loading :SYSTem:PRESet

%     properties(Transient)
    properties
        Ins;
        Set;
    end
    properties
        storage
    end
    methods
        function obj = DSOX3000(addressNUM, xlsName, XlsSheetName)
            if ~exist('addressNUM','var'); addressNUM = 5; end;
            if isnumeric(addressNUM) && round(addressNUM)==addressNUM
                switch addressNUM
                    case  1; address= '0x2A8D::0x1766::MY56311491';
                    case  2; address= '0x0957::0x17A6::MY53160414';
                    case  3; address= '0x0957::0x17A6::MY59121413';
                    case  4; address= '0x0957::0x17A6::MY53160409';
                    case  5; address= '0x2A8D::0x1766::MY54292273';
                    case  6; address= '0x2A8D::0x1766::MY52492279';
                    case  7; address= '0x2A8D::0x1764::MY58262754';
                    case  8; address= '10893::5990::MY56311506';
                    case  9; address= '2391::6054::MY59121408';
                    case 10; address= '0x0957::0x17A6::MY59121311';
                    case 12; address= '0x0699::0x03A6::C046060';
                    case 13; address= '0x0699::0x03A6::C017878';
                    case 14; address= '10893::5988::MY58262754';
                    case 15; address='0x0957::0x17A6::MY53160409';
                    otherwise; error('ERROR! Illegal scope address'); 
%                         obj.Ins =visa('NI', ['USB0::' address '::0::INSTR']);
                end
            elseif ischar(addressNUM); address = addressNUM;
            end
%             obj.Ins = visa('AGILENT', ['USB0::' address '::0::INSTR']);
%             obj.Ins = visa('keysight', ['USB0::' address '::INSTR']);
            obj.Ins = visa('ni', ['USB0::' address '::INSTR']);
%               read_ch = ones(1,4);
            init_flag = false;
            if exist('xlsName','var')
                tic; [NumDat,TxtDat,~] = xlsread(xlsName, XlsSheetName, 'A:E'); toc; 
                obj.Set.Vscale = NumDat(2,1:4); 
                obj.Set.Voffset = NumDat(3,1:4);
                obj.Set.impedance = cell2mat(TxtDat(4,2:5));
                obj.Set.Coupling = TxtDat(5,2:5);
                obj.Set.BWlim = NumDat(6,1:4); 
                obj.Set.N_average = NumDat(7,2); 
                obj.Set.Tscale = NumDat(8,2);
                obj.Set.PreTrigger = NumDat(9,2);
                obj.Set.delay = NumDat(10,2);
                obj.Set.TBref = cell2mat(TxtDat(11,2));
                obj.Set.trig_mode = cell2mat(TxtDat(12,2));
                obj.Set.trig_source = NumDat(13,2);
                obj.Set.trig_slope = NumDat(14,2);
                obj.Set.trig_thresh = NumDat(15,2);
            else 
                obj.Set.Vscale = [0 0 0 0]; % set 0 or negative values to use the current scope setting
                obj.Set.Voffset = [0 0 0 0]; % set values larger than 10 to use the current scope setting
                obj.Set.impedance = ['F' 'F' 'F' 'F']; % 'F' for 50 Ohm, 'M' for 1MegaOhm
                obj.Set.Coupling = {'DC','DC','DC','DC'}; % 'DC' or 'AC'
                obj.Set.BWlim = [0,0,0,0];
                obj.Set.N_average = 1; %32; % set 1 or smaller for non-averaging
                obj.Set.Tscale = 20*1e-3; % set 0 or negative to use current scope settings
                obj.Set.PreTrigger = nan;
                obj.Set.delay = 0;
                obj.Set.TBref = 'CENTer'; % 'CENTer', 'LEFT' or 'RIGHt';
                obj.Set.trig_mode = 'EDGE';
                obj.Set.trig_source = 1; % 'ext' for external;
                obj.Set.trig_slope = 1; % 1 for positive, -1 for negative
                obj.Set.trig_thresh = 0.5;
                obj.Set.Probe = [1,1,1,1];
            end
            
            obj.Ins.InputBufferSize = 2^24;
            if exist('init_flag','var') && init_flag; Init(obj); end;
            
%          fclose(obj.Ins);
        end
     
        function data = IDN(obj) %IDN
            data = query(obj.Ins,'*IDN?');
        end
        
        function Reset(obj)
            fprintf(obj.Ins,'*RST');
        end
        function LoadDefaultSetup(obj)
            fprintf(obj.Ins,':SYSTem:PRESet');
        end
        
        
        
        function Run(obj); obj.isTriggered(); fprintf(obj.Ins, ':RUN'); end;
        function Stop(obj); obj.isTriggered(); fprintf(obj.Ins, ':STOP'); end;
        function Single(obj); obj.isTriggered(); fprintf(obj.Ins, ':SINGle'); end;
        function readyToRead(obj, timeout)
            if exist('timeout', 'var') && islogical(timeout) && timeout; timeout = 11 * obj.getTscale + 0.05; end
            if exist('timeout','var') && ~isempty(timeout) && isnumeric(timeout) && timeout > 0; timeoutclock = tic; end
            while ~obj.isStopped()
                pause(0.05);
                if exist('timeoutclock','var') && toc(timeoutclock) > timeout
                    error('Timeout error. Scope had not stopped');
                end
            end
        end
        function stt = OPC(obj); stt = query(obj.Ins,'*OPC?'); end;
        function stt = OperationComplete(obj); stt = query(obj.Ins,'*OPC?'); end;
        function ClearRegisters(obj); fprintf(obj.Ins,'*CLS?'); end;
        function ForceTrig(obj); fprintf(obj.Ins,':TRIGger:FORCe'); end;
        function stt = Complete(obj); stt = query(obj.Ins,'*OPC?'); end;
        function stt = Status(obj); stt = de2bi(str2double(query(obj.Ins,':OPER:COND?'))); end;
        function stt = isStopped(obj); stt = obj.Status; stt = ~stt(4); end;
        function stt = isTriggered(obj); stt = logical(str2double(query(obj.Ins,'TER?'))); end;
        function [errCode,errString] = readError(obj,disp_flag)
            err = textscan(query(obj.Ins, ':SYSTem:ERRor?'),'%d%q','Delimiter',',');
            errCode = err{1}; errString = err{2}{1};
            if strcmp(errString(1:8),'No error')==1; errCode = false; errString = [];
            elseif exist('disp_flag','var') && disp_flag; disp(errString);
            end
        end
        
        function setVscale(obj, ch, Vscale)
            if length(Vscale) == 1; Vscale = Vscale * ones(1,4); end
            if length(ch) == 4 && all(ch == 1 | ch == 0); ch = find(ch); end
            obj.Set.Vscale(ch) = Vscale(ch);
            for nn = ch
                if Vscale(nn) > 0
                    fprintf(obj.Ins, sprintf(':CHANnel%u:SCALe %g',...
                        nn, obj.Set.Vscale(nn)));
                end
            end
        end
        
        function setVoffset(obj,ch,Voffset)
            if length(Voffset)==1; Voffset = Voffset*ones(1,4); end;
            if length(ch)==4 && all(ch==1 | ch==0); ch = find(ch); end;
            obj.Set.Voffset(ch) = Voffset(ch);
            for nn = ch
                fprintf(obj.Ins, sprintf(':CHANnel%u:OFFSet %g',...
                    nn,obj.Set.Voffset(nn)));
            end
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
        
        function autoVrange(obj,ch,maxTscale,max_tries,force_trig)
            % Tries to optimize VScale to current data. TODO: TEST!!!
            if ~exist('ch','var') || ~isnumeric(ch) || length(ch)~=1 
                error('Please define a single channel for autoscale operation')
            end
            if ~exist('maxTscale','var') || ~isnumeric(maxTscale) || isempty(maxTscale); maxTscale = 10e-3; end
            if ~exist('max_tries','var') || ~isnumeric(max_tries) || isempty(max_tries); max_tries = 10; end
            if ~exist('force_trig','var'); force_trig = 0; end
                                  
            currTscale = obj.getTscale;
            changed_Tscale = 0;
            if currTscale > maxTscale
                changed_Tscale = 1;
                currTdelay = obj.getTdelay;
                obj.setTscale(maxTscale)
                obj.setTdelay(0)
            end
            
            obj.Single;
            if force_trig; obj.ForceTrig; end
            obj.readyToRead;
            [~,v,~] = obj.Read(ch);
            inflag = obj.sigInbounds(v,ch);
            this_try = 1;
            maxVscale = 5;
            this_Vscale = obj.getVscale(ch);
            while ~inflag && (this_try <= max_tries) && (this_Vscale<maxVscale)
                obj.setVscale(ch,this_Vscale*2);
                obj.Single;
                if force_trig; obj.ForceTrig; end
                obj.readyToRead(2);
                [~,v,~] = obj.Read(ch);
                inflag = obj.sigInbounds(v,ch);
                this_try = this_try+1;
                this_Vscale = obj.getVscale(ch);           
                
            end
            if ~inflag
                warning('Failed to do autoVrange, keeping old Vrange')
            else
                obj.setVrange4signal(ch,v,1/1.5)
            end
            if changed_Tscale
                obj.setTscale(currTscale)
                obj.setTdelay(currTdelay)
            end
                
        end
        
        
        function AverageAQN(obj,N_average)
            if N_average>1 && round(N_average)==N_average
                obj.Set.N_average = N_average;
                fprintf(obj.Ins, sprintf(':ACQuire:TYPE %s', 'AVERage')); %% 'NORMal' or 'AVERage'
                fprintf(obj.Ins, sprintf(':ACQuire:COUNt %d', obj.Set.N_average));
%             elseif N_average==1
%                  fprintf(obj.Ins, sprintf(':ACQuire:TYPE %s', 'NORMal')); obj.Set.N_average = 1;
            else obj.Set.N_average = 1; obj.highres; end;
        end
        
        function setChImpedance(obj,ch,impedance)
            if length(impedance)==1; impedance = repmat(impedance,1,4); end;
            if length(ch)==4 && all(ch==1 | ch==0); ch = find(ch); end;
            obj.Set.impedance(ch) = impedance(ch);
            for nn = ch
                if obj.Set.impedance(nn)=='F'
                    fprintf(obj.Ins, sprintf(':CHANnel%d:IMPedance %s',nn,...
                        'FIFTy')); %'FIFTy' or 'ONEMeg'
                elseif obj.Set.impedance(nn)=='M'
                    fprintf(obj.Ins, sprintf(':CHANnel%d:IMPedance %s',nn,...
                        'ONEMeg')); %'FIFTy' or 'ONEMeg'
                end
            end
        end
        
        function setChBW(obj,ch,BW)
%             if iscell(BW);
%                 for Ind = 1:length(BW); temp(Ind,:) = BW{Ind}; end; BW = temp; %#ok<AGROW>
%             end
            if length(BW)==1 && length(ch)>1; BW = repmat(BW,4,1); end;
            if size(BW,1)==1 && size(BW,2)>1; BW = BW.'; end;
            if length(ch)==4 && all(ch==1 | ch==0); ch = find(ch); end;
            obj.Set.BWlim(ch) = BW(ch);
            for nn = ch
                fprintf(obj.Ins, sprintf(':CHANnel%d:BWLimit %d',nn,obj.Set.BWlim(nn)));
            end
        end
        
        function setChCoupling(obj,ch,coupling)
            if iscell(coupling);
                for Ind = 1:length(coupling); temp(Ind,:) = coupling{Ind}; end; %#ok<AGROW>
                coupling = temp;
            end
            if size(coupling,1)==1 && length(ch)>1
                coupling = repmat(coupling,4,1);
            end
            for nn = ch
                if length(ch)==1
                    obj.Set.Coupling{nn} = coupling;
                    fprintf(obj.Ins, sprintf(':CHANnel%d:COUPling %s',nn,coupling));
                else
                    obj.Set.Coupling{nn} = coupling(nn,:);
                    fprintf(obj.Ins, sprintf(':CHANnel%d:COUPling %s',nn,coupling(nn,:)));
                end
            end
        end
        
        function setTscale(obj,Tscale)
            obj.Set.Tscale = Tscale;
            if Tscale>0
                fprintf(obj.Ins, sprintf(':TIMebase:SCALe %g', obj.Set.Tscale));
            end
        end
        function T = getTscale(obj); T = str2double(query(obj.Ins, ':TIMebase:SCALe?')); end

        function setTdelay(obj,delay)
            obj.Set.delay = delay;
            fprintf(obj.Ins, sprintf(':TIMebase:POSition %g', obj.Set.delay));
        end
        function T=getTdelay(obj); T=str2double(query(obj.Ins,':TIMebase:POSition?')); end;

        function setTref(obj,ref)
            if ~exist('ref','var')
                ref = 'LEFT'; % 'CENTer'; % 'RIGHt';
%                 ref = questdlg('Choose timebase reference position',...
%                     'timebase reference','LEFT','CENTer','RIGHt','LEFT');
            end
            obj.Set.TBref = ref;
            fprintf(obj.Ins, sprintf([':TIMebase:REFerence ' obj.Set.TBref]));
        end

        function TrigThresh(obj, thresh)
            % setting the trigger threshhold in volts
            obj.Set.trig_thresh = thresh;
            fprintf(obj.Ins, sprintf(':TRIGger:EDGE:LEVel %g', obj.Set.trig_thresh));
        end
        
        function SrcOut = TrigSource(obj, src)
            if ~exist('src','var'); SrcOut = query(obj.Ins,':TRIGger:EDGE:SOURce?'); return; end
            SrcOut = src; obj.Set.trig_source = src;
            if isnumeric(src); src = sprintf('CHANnel%d',src);
            elseif strcmpi(src(1:3),'ext'); src = 'EXTernal';
            elseif strcmpi(src(1:3),'lin'); src = 'LINE';
            end
            fprintf(obj.Ins, sprintf(':TRIGger:EDGE:SOURce %s', src));
        end
        
        function TrigSlope(obj,slope)
            if isnumeric(slope)
                if slope>0
                    fprintf(obj.Ins, sprintf(':TRIGger:EDGE:SLOPe %s', 'POSitive'));
                elseif slope<0
                    fprintf(obj.Ins, sprintf(':TRIGger:EDGE:SLOPe %s', 'NEGative'));
                end
            elseif ischar(slope)
                    fprintf(obj.Ins, sprintf(':TRIGger:EDGE:SLOPe %s', slope));
            else error('ERROR! Incorrect SLOPE trigger definition');
            end
        end
        
        function TrigMode(obj,swe)
            fprintf(obj.Ins, sprintf(':TRIG:SWE %s', swe));
        end
        
        function [t,v,sat_flag] = Read(obj,ReadChannels)
%             if ~obj.isTriggered; t = []; v = []; return; end;
%             if ~obj.isStopped; obj.Stop; end;
            obj.OPC();
%             if ~logical(str2double(query(obj.Ins,':WAV:POIN?'))); t = []; v = []; return; end;
            if length(ReadChannels)==4 && all(ReadChannels==1 | ReadChannels==0)
                temp = 1:4; ReadChannels = temp(find(ReadChannels)); %#ok<FNDSB>
            end
            for ch = ReadChannels
                fprintf(obj.Ins, [':WAV:SOUR CHAN', num2str(ch)]);
        		fprintf(obj.Ins, ':WAV:BYT LSBF');
                fprintf(obj.Ins, ':WAV:FORM WORD');
%                 fprintf(obj.Ins, ':WAV:POIN:MODE MAX');
%                 fprintf(obj.Ins, ':WAV:POIN 62500');
		
        		WFMPRE = strsplit(query(obj.Ins,':WAV:PRE?'),','); obj.OPC;
                NumPoints = str2double(WFMPRE{3});
                XIncr = str2double(WFMPRE{5});
                XOrig = str2double(WFMPRE{6});
                XRef = str2double(WFMPRE{7});
                YIncr = str2double(WFMPRE{8});
                YOrig = str2double(WFMPRE{9});
                YRef = str2double(WFMPRE{10});
                
                if ch==ReadChannels(1)
        			v = zeros(length(ReadChannels),NumPoints);
                    t = ((1:NumPoints) - XRef)*XIncr + XOrig;
                    sat_flag = false(length(ReadChannels),1);
                end
                
                fprintf(obj.Ins, ':WAV:DATA?');
                fread(obj.Ins,10);
                V = fread(obj.Ins, NumPoints, 'uint16');
                if ~strcmp(char(fread(obj.Ins,1)),sprintf('\n')); error('ERROR! number of points in waveform is not as expected'); end; %#ok<FREAD>
                sat_flag(ReadChannels==ch) = any(diff(V==1)==1) || ...
                    all(V==1)==1 || any(diff(V==1)==65535) || all(V==1)==65535;
				V = reshape(V, [1 NumPoints]);
        		v(find(ReadChannels==ch),:) = (V - YRef)*YIncr + YOrig; %#ok<FNDSB>
            end
        end
        
        function offset = getVoffset(obj,ch)
            if ~exist('ch','var'); ch = 1:4; end;
            offset = nan(4,1);
            for nn = ch
                fprintf(obj.Ins,sprintf(':CHAN%d:OFFS?',nn));
                if length(ch)==1; offset = str2double(fscanf(obj.Ins));
                else offset(nn) = str2double(fscanf(obj.Ins));
                end
            end
        end
        
        function scale = getVscale(obj, ch)
            % get the voltage per devision scale of channel
            if ~exist('ch', 'var'); ch = 1:4; end
            scale = nan(4, 1);
            for nn = ch
                fprintf(obj.Ins,sprintf(':CHAN%d:SCAL?',nn));
                if length(ch)==1; scale = str2double(fscanf(obj.Ins));
                else
                    scale(nn) = str2double(fscanf(obj.Ins));
                end
            end
        end
        
        function R = getVrange(obj, ch)
            % get range of signal in ch
            if ~exist('ch', 'var'); ch = 1:4; end
            V0 = getVoffset(obj, ch); dV = getVscale(obj, ch);
            R = [V0(:) - 5 * dV(:), V0(:) + 5 * dV(:)];
        end
        
        function Init(obj)
            fprintf(obj.Ins,':TRIGger:ZONE:STATe 0');
            obj.ClearRegisters;
            init_cahnnels(obj,1:4);
            init_horizontal(obj);
            init_trig(obj);
            AverageAQN(obj,obj.Set.N_average);
            if obj.Set.N_average<2; fprintf(obj.Ins, sprintf(':ACQuire:TYPE %s', 'HRESolution')); end;
        end
        function init_cahnnels(obj,ch)
            if ~exist('ch','var'); ch = 1:4; end;
            if length(ch)>1; setChCoupling(obj,ch,{obj.Set.Coupling{1:4}});
            else setChCoupling(obj,ch,obj.Set.Coupling{ch});
            end
            setChBW(obj,ch,obj.Set.BWlim(1:4));
            setChImpedance(obj,ch,obj.Set.impedance(1:4));
            setVscale(obj,ch,obj.Set.Vscale(1:4));
            setVoffset(obj,ch,obj.Set.Voffset(1:4));
            Probe(obj,ch,obj.Set.Probe(1:4));
%             if isfield(params,'BW'); setChBW(obj,ch,params.BW); end;            
%             if isfield(params,'impedance'); setChImpedance(obj,ch,params.impedance); end;
%             if isfield(params,'Vscale'); setVscale(obj,ch,params.Vscale); end;
%             if isfield(params,'Voffset'); setVoffset(obj,ch,params.Voffset); end;
        end
        function init_horizontal(obj)
            fprintf(obj.Ins, sprintf(':TIMebase:REFerence %s', obj.Set.TBref));
            setTscale(obj,obj.Set.Tscale);
            setTdelay(obj,obj.Set.delay);
%             if isfield(params,'Tscale'); setTscale(obj,params.Tscale); end;
%             if isfield(params,'delay'); setTdelay(obj,params.delay); end;
        end
        function init_trig(obj)
%             if ~exist('params','var'); params = []; end;
            
              fprintf(obj.Ins, sprintf(':TRIGger:MODE %s', obj.Set.trig_mode));
              TrigSource(obj,obj.Set.trig_source);
              TrigThresh(obj,obj.Set.trig_thresh);
              
              if obj.Set.trig_slope>0
                  fprintf(obj.Ins, sprintf(':TRIGger:EDGE:SLOPe %s', 'POSitive'));
              elseif obj.Set.trig_slope<0
                  fprintf(obj.Ins, sprintf(':TRIGger:EDGE:SLOPe %s', 'NEGative'));
              end

%             if isfield(params,'trig_mode')
%                 fprintf(obj.Ins, sprintf(':TRIGger:MODE %s', params.trig_mode));
%             end
            
%             if isfield(params,'trig_source'); TrigSource(obj,params.trig_source); end;
%             if isfield(params,'trig_thresh'); TrigThresh(obj,params.trig_thresh); end;
            
%             if isfield(params,'trig_slope')
%                 if params.trig_slope>0
%                     fprintf(obj.Ins, sprintf(':TRIGger:EDGE:SLOPe %s', 'POSitive'));
%                 elseif params.trig_slope<00
%                     fprintf(obj.Ins, sprintf(':TRIGger:EDGE:SLOPe %s', 'NEGative'));
%                 end
%             end
            
        end
        
        function updated = updateSet(obj)
            for ch = 1:4
                updated.Vscale(ch) = str2double(query(obj.Ins,[':CHANnel' num2str(ch) ':SCAL?']));
                updated.Voffset(ch) = str2double(query(obj.Ins,[':CHANnel' num2str(ch) ':OFFS?']));
                switch query(obj.Ins,[':CHANnel' num2str(ch) ':IMP?'])
                    case sprintf('ONEM\n'); updated.impedance(ch) = 'M';
                    case sprintf('FIFT\n'); updated.impedance(ch) = 'F';
                end
                updated.Coupling{ch} = query(obj.Ins,[':CHANnel' num2str(ch) ':COUP?']);
                updated.Coupling{ch} = updated.Coupling{ch}(1:end-1);
                updated.BWlim(ch) = str2double(query(obj.Ins,[':CHANnel' num2str(ch) ':BWLimit?']));
            end
            if strcmpi(query(obj.Ins,':ACQuire:TYPE?'),sprintf('NORM\n'))
                updated.N_average = 1;
            else updated.N_average = str2double(query(obj.Ins,':ACQuire:COUNt?'));
            end
            updated.Tscale = str2double(query(obj.Ins,':TIMebase:SCALe?'));
            updated.delay = str2double(query(obj.Ins,':TIMebase:POSition?'));
            updated.PreTrigger = nan;
            updated.TBref = query(obj.Ins,':TIMebase:REFerence?'); updated.TBref = updated.TBref(1:end-1);
            updated.trig_mode = query(obj.Ins,':TRIGger:MODE?'); updated.trig_mode = updated.trig_mode(1:end-1);
            updated.trig_source = query(obj.Ins,':TRIGger:SOUR?'); updated.trig_source = updated.trig_source(1:end-1);
            switch query(obj.Ins,':TRIGger:EDGE:SLOP?')
                case sprintf('POS\n'); updated.trig_slope = 1;
                case sprintf('NEG\n'); updated.trig_slope = -1;
                otherwise; updated.trig_slope = nan;
            end
            updated.trig_thresh = str2double(query(obj.Ins,':TRIGger:EDGE:LEV?'));
        end
        
        function stt_out = Probe(obj,ch,stt_in)
            stt_out = zeros(size(ch));
            for c = ch
                if exist('stt_in','var')
                    fprintf(obj.Ins,[':CHANnel' num2str(c) ':PROBe ' num2str(stt_in(c==ch))]);
                end
                stt_out(c==ch) = str2double(query(obj.Ins,[':CHANnel' num2str(c) ':PROBe?']));
            end
        end
        
        function setChDisp(obj,ch); fprintf(obj.Ins,[':CHANnel' num2str(ch) ':DISPlay 1']); end
        
        function highres(obj); fprintf(obj.Ins, sprintf(':ACQuire:TYPE %s', 'HRESolution')); end;
        
        function Z = getChImpedance(obj,ch)
            Z = zeros(size(ch));
            for c=ch
                answer = query(obj.Ins,sprintf(':CHANnel%d:IMPedance?',c));
                switch lower(answer(1)); case 'f'; Z(ch==c) = 50; case 'o'; Z(ch==c) = 1e6; end
            end
        end
        
        function Z = getChCoupling(obj,ch)
            Z = cell(size(ch));
            for c=ch
                Z{ch==c} = query(obj.Ins,sprintf(':CHANnel%d:COUPling?',c));
                Z{ch==c} = Z{ch==c}(1:2);
            end
        end
        
        function outtext =  get_channel_setting(obj, ch)
            fprintf(obj.Ins,[':CHANnel' num2str(ch) ':SCALe?']);vscale = fscanf(obj.Ins);
            fprintf(obj.Ins,[':CHANnel' num2str(ch) ':OFFSet?']);voffset = fscanf(obj.Ins);
            fprintf(obj.Ins,[':CHANnel' num2str(ch) ':COUPling?']);coupling = fscanf(obj.Ins);
            fprintf(obj.Ins,[':CHANnel' num2str(ch) ':IMPedance?']);impedance = fscanf(obj.Ins);
            
            outtext = ['CH: ' num2str(ch) ' ' 'Vscale: ' num2str(vscale) ];
            outtext = [outtext ' ' 'Voffset: ' num2str(voffset)];
            outtext = [outtext ' ' 'Coupling: ' coupling];
            outtext = [outtext ' ' 'Impedance: ' impedance];
        end
        
        function [t,v,SatFlag] = ReadWFM(obj,ReadChannels)
            obj.OPC();
            if length(ReadChannels)==4 && all(ReadChannels==1 | ReadChannels==0)
                temp = 1:4; ReadChannels = temp(find(ReadChannels)); %#ok<FNDSB>
            end
            SatFlag = false(size(ReadChannels));
            for ch = ReadChannels

                fprintf(obj.Ins, [':WAV:SOUR CHAN', num2str(ch)]);
        		fprintf(obj.Ins, ':WAV:BYT LSBF');
                fprintf(obj.Ins, ':WAV:FORM WORD');
        		WFMPRE = strsplit(query(obj.Ins,':WAV:PRE?'),','); obj.OPC;
                NumPoints = str2double(WFMPRE{3});
                XIncr = str2double(WFMPRE{5});
                XOrig = str2double(WFMPRE{6});
                XRef = str2double(WFMPRE{7});
                YIncr = str2double(WFMPRE{8});
                YOrig = str2double(WFMPRE{9});
                YRef = str2double(WFMPRE{10});
                
                if ch==ReadChannels(1)
        			v = zeros(length(ReadChannels),NumPoints);
                    t = ((1:NumPoints) - XRef)*XIncr + XOrig;
                end
                
                fprintf(obj.Ins, ':WAV:DATA?');
                fread(obj.Ins, 10);
                V = fread(obj.Ins, NumPoints, 'uint16');
				V = reshape(V, [1 NumPoints]);
                SatFlag(ReadChannels==ch) = sum(abs(V)==256)>(0.05*length(V)) ...
                    || sum(abs(V)==65280)>(0.05*length(V));
        		v(find(ReadChannels==ch),:) = (V - YRef)*YIncr + YOrig; %#ok<FNDSB>
            end
        end
        function HighRes(obj); fprintf(obj.Ins,':ACQuire:TYPE HRES'); end
        function Normal(obj); fprintf(obj.Ins,':ACQuire:TYPE NORM'); end

        function setCurrentProbe(obj,ch,ratio)
            if ~exist('ratio','var'); ratio = 1; end
            if ischar(ratio) && ( strcmpi(ratio,'Rogovski') || strcmpi(ratio,'Peerson') ); ratio = 10; end
            if ischar(ratio) && ( strcmpi(ratio,'CT2') ); ratio = 1; end
            if ischar(ratio) && ( strcmpi(ratio,'CT1') ); ratio = 1/5; end
            if length(ratio)==1 && length(ch)>1; ratio = taio*ones(size(ch)); end
            for c = ch
                fprintf(obj.Ins, sprintf(':CHANnel%d:UNITs %s',c,'AMPere'));
                fprintf(obj.Ins, sprintf(':CHANnel%d:PROBe X%d',c,ratio(c==ch)));
                obj.setChImpedance(c,'F');
            end
        end

        
        function flag = sigAdjusted(obj,S,ch,rng,tol)
            if ~exist('rng','var') || isempty(rng) || ~(rng<=1 && rng>=0); rng = 0.6; end
            if ~exist('tol','var') || isempty(tol) || ~(tol<=1 && tol>=0); tol = 0.1; end
            R = obj.getVrange(ch);
            flag = abs(max(S)-((0.5+rng/2)*max(R)+(0.5-rng/2)*min(R)))<(tol*diff(R)) && ...
                   abs(min(S)-((0.5-rng/2)*max(R)+(0.5+rng/2)*min(R)))<(tol*diff(R));
            % Gil's addition 2.7
            inboundsflag = max(S) < R(2) && min(S) > R(1);
            flag = flag && inboundsflag;
        end
        function inboundsflag = sigInbounds(obj, sig, ch)
            % checks if signal contained in scope bounds
            scope_bounds = obj.getVrange(ch);
            inboundsflag = (max(sig) < scope_bounds * [0.02; 0.98] ) && (min(sig) > scope_bounds * [0.98; 0.02]);
        end
        function minscaleflag = sigMinScale(obj,S,ch,minSscale)
            R = obj.getVrange(ch);
            minscaleflag = ((max(S)-min(S))/diff(R))>minSscale;
        end
        
        function Nout = Npoints(obj,N)
            if ~exist('N','var') || isempty(N) || ~isnumeric(N) || isnan(N)
                Nout = num2str(query(obj.Ins,':WAV:POIN?')); return;
            elseif isinf(N)
                fprintf(obj.Ins,':WAV:POIN MAX');
            elseif N<=0 || N~=round(N)
                error('ERROR! Can only ask for a round positive number of points');
            else; fprintf(obj.Ins,[':WAV:POIN ' num2str(N)]);
            end
            obj.OPC;
            Nout = num2str(query(obj.Ins,':WAV:POIN?'));
            obj.OPC;
            obj.readError(true);
        end
   end
end