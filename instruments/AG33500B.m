classdef AG33500B
    % Updates:
    % 10/03/19 - Roy: Added VLim function (limits the voltage output) - line 254
    % 16/04/19 - Roy: Added Output function (turns\reads output state) - line 65
    % 15/06/19 - Roy: Added SingleTrapz function + cosmetic changes - lines 37, 192
    % 12/07/19 - Roy: BurstSin changes: Fixed trigger bug + added delay - lines 112,122,126
    % 26/07/19 - Roy: BurstON & Burst functions added + lines 290-304
    % 30/08/19 - Roy: Added Delay function + lines 333-347
    % 12/09/19 - Avraham Berrebi: Added PulseW for pulse with specific rise and fall time
    %                             Added BurstON1 - for 1 cycle, BurstON_INF - for infinity cycles
    % 26/01/20 - Roy: Added Abort function - lines 404
    % 02/02/20 - Roy: Added setTrig function - line 408
    % 13/02/20 - Roy: Added "all_waveform_data" flag to Output function - line 86
    % 13/02/20 - Roy: Added readError function - line 51
    % 05/03/20 - Roy: Added SinglePulseV2 function which replaces singlePulse - line 261
    % 05/03/20 - Roy: Added BurstSinV2 function which replaces BurstSin - line 261
    % 11/03/20 - Roy: Added "Burst-off" to DC - line 133
    % 08/07/20 - Gil: Sin function did not work, pasted code from old version
    % and fixed it.
    % 05/08/20 - Itay: Added SetTrigger, which is the same as SetTrig
    % 22/09/20 - Alon: Added SweepConfig, SweepON, SweepOFF
    % 26/07/21 - Elkabetz: Added function getChannelStateNum which returns the state in nums
    % 31/03/22 - CF: Added output termination (impedance) function
    
    %    properties(Transient)
    properties
        Ins
        DefARBname
    end
    methods
        function obj = AG33500B(serial,connection_type)
            if ~exist('serial','var'); serial = '02252'; end
            if ~exist('connection_type','var'); connection_type = 'tcpip'; end
            switch lower(connection_type)
                case 'tcpip'
                    if isnumeric(serial); serial = num2str(serial); end
                    obj.Ins = instrfind('Type','visa-tcpip','name',['TCPIP0::' serial '::inst0::INSTR']);
                    if isempty(obj.Ins)
                        obj.Ins = visa('NI',['TCPIP0::' serial '::inst0::INSTR']);
                    else
                        fclose(obj.Ins);
                        obj.Ins = obj.Ins(1);
                    end
                case 'usb'
                    if strcmp(serial(1:3),'USB')
                        obj.Ins = visa('NI', serial);
                    else
                        obj.Ins = visa('NI', ['USB0::0x0957::' serial '::0::INSTR']);
                    end
            end
            set(obj.Ins,'InputBufferSize',2^18);
            set(obj.Ins,'OutputBufferSize',2^18);
            obj.Ins.Timeout = 2;
            
            fopen(obj.Ins);
            obj.Polarity([1,2],+1);
            
            
            obj.DefARBname = 'AWFM34';
            if ~strcmp(serial(1:3),'USB'); fclose(obj.Ins); end
        end
        
        function data = IDN(obj); data = query(obj.Ins,'*IDN?'); end %IDN
        
        function [errCode,errString] = readError(obj,disp_flag)
            err = textscan(query(obj.Ins, ':SYSTem:ERRor?'),'%d%q','Delimiter',',');
            errCode = err{1}; errString = err{2}{1};
            if strcmp(errString(1:8),'No error')==1; errCode = false; errString = [];
            elseif exist('disp_flag','var') && disp_flag; disp(errString);
            end
        end
        
        function Reset(obj); fprintf(obj.Ins, '*RST'); end
        
        function Trig(obj); fprintf(obj.Ins, '*TRG'); end
        
        function out = OPC(obj); out = str2double(query(obj.Ins,'*OPC?')); end
        
        function beep(obj); fprintf(obj.Ins,'SYSTem:BEEPer:IMMediate'); end
        
        function OutputON(obj,ch)
            if ~exist('ch','var'); ch = 1:2; end
            if any(ch==1); fprintf(obj.Ins, 'OUTP1 1'); end
            if any(ch==2); fprintf(obj.Ins, 'OUTP2 1'); end
        end
        
        function OutputOFF(obj,ch)
            if ~exist('ch','var'); ch = 1:2; end
            if any(ch==1); fprintf(obj.Ins, 'OUTP1 0'); end
            if any(ch==2); fprintf(obj.Ins, 'OUTP2 0'); end
        end
        function Termination(obj,ch,imp)
            if ~exist('ch','var'); ch = 1:2; end
            if any(ch==1); fprintf(obj.Ins, ['OUTP1:LOAD ' num2str(imp)]); end
            if any(ch==2); fprintf(obj.Ins, ['OUTP2:LOAD ' num2str(imp)]); end
        end
        
        function st = OutputState(obj, ch)
            if ~exist('ch', 'var'); ch = 1:2; end
            if any(ch == 1)
                st1 = str2double(query(obj.Ins, 'OUTP1?'));
            else; st1 = [];
            end
            if any(ch == 2)
                st2 = str2double(query(obj.Ins, 'OUTP2?'));
            else; st2 = [];
            end
            st = [st1, st2];
        end
        
        function [stt, wvfm] = Output(obj,ch,stt,all_waveform_data)
            if ~exist('ch', 'var') || isempty(ch) || ~isnumeric(ch) || (any(ch ~= 1) && any(ch ~= 2)); ch = 1:2; end
            if ~exist('stt', 'var') || isempty(stt) || ~isnumeric(stt) || ~islogical(stt)
                for I = ch
                    stt(I == ch) = str2double(query(obj.Ins,['OUTP' num2str(I) '?']));
                end
            else
                for I = ch
                    fprintf(obj.Ins,['OUTP' num2str(I) ' ' num2str(stt(I == ch))]);
                end
            end
            if exist('all_waveform_data', 'var') && islogical(all_waveform_data) && all_waveform_data
                for I = ch
                    temp = textscan(query(obj.Ins,['SOUR' num2str(I) ':APPL?']),'%s %f,%f,%f');
                    wvfm{I == ch}.name = temp{1}{1}(2:end); for J = 2:length(temp); wvfm{I==ch}.p(J-1)= temp{J}(1); end
                end
            else; wvfm = [];
            end
        end

        function data = getPhase(obj,ch)
            data = str2double(query(obj.Ins, ['SOURce' num2str(ch) ':PHAS?']));
        end
        
        function Sin(obj, ch, freq, Amp, phase, offset)
            ch = num2str(ch); freq = num2str(freq); Amp = num2str(Amp);
            if ~exist('phase','var') || isempty(phase) || isnan(phase); phase = 0; end; 
            if ~exist('offset','var') || isempty(offset) || isnan(offset); offset = 0; end; offset = num2str(offset);
            
            fprintf(obj.Ins,['SOURce' ch ':FUNC SIN']);
            fprintf(obj.Ins,['SOURce' ch ':FREQ ' freq]);
            fprintf(obj.Ins,['SOURce' ch ':VOLT ' Amp]);
            fprintf(obj.Ins,['SOURce' ch ':VOLT:OFFS ' offset]);

            if obj.getPhase(str2double(ch)) ~= phase
                fprintf(obj.Ins,['SOURce' ch ':PHAS ' num2str(phase)]);
            end
        end
        
        function Square(obj, ch, freq, Amp, offset, phase)
            ch = num2str(ch); freq = num2str(freq);
            if ~exist('offset','var'); offset = Amp; end; offset = num2str(offset);
            if ~exist('phase','var'); phase = 0; end; 
            Amp = num2str(2*Amp);
            
            fprintf(obj.Ins,['SOURce' ch ':FUNC SQUare']);
            fprintf(obj.Ins,['SOURce' ch ':FREQ ' freq]);
            fprintf(obj.Ins,['SOURce' ch ':VOLT ' Amp]);
            fprintf(obj.Ins,['SOURce' ch ':VOLT:OFFS ' offset]);
            if obj.getPhase(str2double(ch)) ~= phase
                fprintf(obj.Ins,['SOURce' ch ':PHAS ' num2str(phase)]);
            end
            
        end
        
        function DC(obj, ch, value)
            ch = num2str(ch); value = num2str(value);
            if str2double(query(obj.Ins,['SOURce' ch ':BURSt:STATe?'])); fprintf(obj.Ins,['SOURce' ch ':BURSt:STATe OFF']); end
            fprintf(obj.Ins,['SOURce' ch ':FUNC DC']);
            fprintf(obj.Ins,['SOURce' ch ':VOLT:OFFS ' value]);
        end
        
        % THIS FUNCITION IS OUTDATED! USE BurstSinV2 INSTEAD.
        function BurstSin(obj,ch,Ncyc,trig,freq,amp,offset,delay)
            % THIS FUNCITION IS OUTDATED! USE BurstSinV2 INSTEAD.
            if ~exist('ch','var') || isempty(ch); ch = 1; end
            if ~exist('offset','var') || isempty(offset); offset = 0; end
            if ~exist('trig','var') || isempty(trig); trig = 'BUS'; end
            if ~exist('delay','var') || isempty(delay); delay = 0; end
            src = [':SOURce' num2str(ch)];
            
            fprintf(obj.Ins,[src ':APPLy:SINusoid ' ...
                num2str(freq) ',' num2str(amp) ',' num2str(offset) ]);
            fprintf(obj.Ins, [src ':BURSt:STATe ON']);
            ch = num2str(ch);
            fprintf(obj.Ins,['TRIGger' num2str(ch) ':SOURce ' trig]);
            fprintf(obj.Ins,[src ':BURSt:NCYC ' num2str(Ncyc)]);
            fprintf(obj.Ins,[src ':BURSt:MODE TRIGgered']);
            fprintf(obj.Ins,['TRIGger' num2str(ch) ':DELay ' num2str(delay)]);
            
            fprintf(obj.Ins,'OUTPut:TRIGger ON');
            fprintf(obj.Ins,['OUTPut:TRIGger:SOURce CH' num2str(ch)]);
            %         fprintf(obj.Ins,['TRIGger' num2str(ch) ':SOURce BUS']);
            fprintf(obj.Ins,['OUTPut:SYNC:SOURce CH' num2str(ch)]);
        end
        
        function Pulse(obj,ch,freq,V,dt,offset,phase)
            if ~exist('ch','var') || isempty(ch) || isnan(ch); ch = 1; end
            if ~exist('freq','var') || isempty(freq) || isnan(freq); freq = 1/(4*dt); end
            if ~exist('offset','var') || isempty(offset) || isnan(offset); offset = 0; end
            if ~exist('phase','var') || isempty(phase) || isnan(phase); phase = nan; end
            src = [':SOURce' num2str(ch)];
            
            fprintf(obj.Ins,[src ':APPLy:PULSe ' num2str(freq) ',' num2str(V) ',' num2str(offset+V/2)]);
            fprintf(obj.Ins,[src ':FUNCtion:PULSe:WIDTh ' num2str(dt)]);
            if ~isnan(phase)
                fprintf(obj.Ins,[src ':PHASe:SYNC']);
                fprintf(obj.Ins,[src ':PHASe ' num2str(phase)]);
            end
        end
        
        function PulseW(obj,ch,freq,amp,delay,transL,transT,width)
            src = [':SOURce' num2str(ch)];
            fprintf(obj.Ins,[src ':APPLy:PULSe ' num2str(freq) ',' num2str(amp) ',' num2str(amp/2)]);
            fprintf(obj.Ins,[src ':FUNCtion:PULSe:WIDTh ' num2str(width)]);
            fprintf(obj.Ins,['TRIGger',num2str(ch),':SOURce ','EXTernal']);
            fprintf(obj.Ins,['SOUR',num2str(ch),':BURSt:STATe ',num2str(1)']); %always triggered
            fprintf(obj.Ins,['SOUR',num2str(ch),':FUNC:PULS:TRAN:LEAD ',num2str(transL)]);
            fprintf(obj.Ins,['SOUR',num2str(ch),':FUNC:PULS:TRAN:TRA +',num2str(transT)]);
            fprintf(obj.Ins,['TRIGger',num2str(ch),':DELay ',num2str(delay)]);
            fprintf(obj.Ins,['OUTP',num2str(ch),' ',num2str(1)]); % aways on
        end
        
        function BurstPulse(obj,ch,freq,V,dt,offset,delay) % dt is in Sec (length of each pulse)
            %         if ~exist('phase','var') || isempty(phase) || isnan(phase); phase = nan; end
            src = [':SOURce' num2str(ch)];
            
            fprintf(obj.Ins,[src ':APPLy:PULSe ' num2str(freq) ',' num2str(V) ',' num2str(offset+V/2)]);
            fprintf(obj.Ins,[src ':FUNCtion:PULSe:WIDTh ' num2str(dt)]);
            %         if ~isnan(phase)
            %             fprintf(obj.Ins,[src ':PHASe:SYNC']);
            %             fprintf(obj.Ins,[src ':PHASe ' num2str(phase)]);
            %         end
            fprintf(obj.Ins,['TRIGger',num2str(ch),':SOURce EXTernal']);
            fprintf(obj.Ins,[src ':BURSt:STATe ON']);
            fprintf(obj.Ins,[src ':BURSt:NCYC 1']);
            fprintf(obj.Ins,[src ':BURSt:MODE TRIGgered']);
            
            %         fprintf(obj.Ins,'OUTPut:TRIGger ON');
            %         fprintf(obj.Ins,['OUTPut:TRIGger:SOURce CH' num2str(ch)]);
            fprintf(obj.Ins,['OUTPut:SYNC:SOURce CH' num2str(ch)]);
            
            if exist('delay','var')
                fprintf(obj.Ins,['TRIGger',num2str(ch),':DELay ',num2str(delay)]);
            end
        end
        
        function Ramp(obj,ch,freq,V,offset,phase)
            if ~exist('ch','var') || isempty(ch) || isnan(ch); ch = 1; end
            if ~exist('offset','var') || isempty(offset) || isnan(offset); offset = 0; end
            if ~exist('phase','var') || isempty(phase) || isnan(phase); phase = nan; end
            src = [':SOURce' num2str(ch)];
            fprintf(obj.Ins,[src ':APPLy:RAMP ' num2str(freq) ',' num2str(V) ',' num2str(offset)]);
            if ~isnan(phase)
                fprintf(obj.Ins,[src ':PHASe:SYNC']);
                fprintf(obj.Ins,[src ':PHASe ' num2str(phase)]);
            end
        end
        
        function Triangle(obj,ch,freq,V,offset,phase)
            if ~exist('ch','var') || isempty(ch) || isnan(ch); ch = 1; end
            if ~exist('offset','var') || isempty(offset) || isnan(offset); offset = 0; end
            if ~exist('phase','var') || isempty(phase) || isnan(phase); phase = nan; end
            src = [':SOURce' num2str(ch)];
            fprintf(obj.Ins,[src ':APPLy:TRIangle ' num2str(freq) ',' num2str(V) ',' num2str(offset)]);
            if ~isnan(phase)
                fprintf(obj.Ins,[src ':PHASe:SYNC']);
                fprintf(obj.Ins,[src ':PHASe ' num2str(phase)]);
            end
        end
        
        % THIS FUNCITION IS OUTDATED! USE SinglePulseV2 INSTEAD.
        function SinglePulse(obj,ch,V,dt,TrigSource,offset,dt0)
            % THIS FUNCITION IS OUTDATED! USE SinglePulse_v2 INSTEAD.
            if ~exist('ch','var') || isempty(ch) || isnan(ch); ch = 1; end
            if ~exist('TrigSource','var'); TrigSource = 'BUS'; end
            if ~exist('offset','var'); offset = 0; end
            if ~exist('dt0','var') || isempty(dt0) || isnan(dt0); dt0 = 0; end
            src = [':SOURce' num2str(ch)];
            
            if V<offset && obj.Polarity(ch)>0; obj.Polarity(ch,-1); end
            if V>offset && obj.Polarity(ch)<0; obj.Polarity(ch,+1); end
            fprintf(obj.Ins,[src ':APPLy:PULSe ' num2str(1/(2*dt)) ',' num2str(abs(V-offset)) ',' num2str((V+offset)/2)]);
            fprintf(obj.Ins,[src ':FUNCtion:PULSe:WIDTh ' num2str(dt)]);
            fprintf(obj.Ins,[src ':BURSt:STATe ON']);
            fprintf(obj.Ins,['TRIGger' num2str(ch) ':SOURce BUS']);
            fprintf(obj.Ins,['TRIGger' num2str(ch) ':DELay ' num2str(dt0)]);
            fprintf(obj.Ins,[src ':BURSt:NCYC 1']);
            fprintf(obj.Ins,[src ':BURSt:MODE TRIGgered']);
            OutputON(obj,ch);
            
            fprintf(obj.Ins,'OUTPut:TRIGger ON');
            fprintf(obj.Ins,['OUTPut:TRIGger:SOURce CH' num2str(ch)]);
            fprintf(obj.Ins,['TRIGger' num2str(ch) ':SOURce ' TrigSource]);
            fprintf(obj.Ins,['OUTPut:SYNC:SOURce CH' num2str(ch)]);
            
            if strcmpi(TrigSource,'BUS'); Trig(obj); end
        end
        
        function SinglePulseV2(obj,ch,V,dt,TrigSource,offset,dt0)
            if ~exist('ch','var') || isempty(ch) || isnan(ch); ch = 1; end
            if ~exist('TrigSource','var'); TrigSource = 'BUS'; end
            if ~exist('offset','var'); offset = 0; end
            if ~exist('dt0','var') || isempty(dt0) || isnan(dt0); dt0 = 0; end
            src = [':SOURce' num2str(ch)];
            
            if V<offset && obj.Polarity(ch)>0; obj.Polarity(ch,-1); end
            if V>offset && obj.Polarity(ch)<0; obj.Polarity(ch,+1); end
            [~,wvfm] = obj.Output(ch,nan,true); if ~strcmp(wvfm{1}.name,'PULS'); fprintf(obj.Ins,[src ':FUNCtion PULSe']); end
            
            fprintf(obj.Ins,[src sprintf(':FREQuency %g', 1/(2*dt))]);
            fprintf(obj.Ins,[src sprintf(':VOLTage %g', abs(V-offset))]);
            fprintf(obj.Ins,[src sprintf(':VOLTage:OFFSet %g', (V+offset)/2)]);
            
            %         fprintf(obj.Ins,['TRIGger' num2str(ch) ':SOURce ' TrigSource]);
            fprintf(obj.Ins,[src ':BURSt:STATe ON']);
            fprintf(obj.Ins,[src ':FUNCtion:PULSe:WIDTh ' num2str(dt)]);
            fprintf(obj.Ins,['TRIGger' num2str(ch) ':DELay ' num2str(dt0)]);
            fprintf(obj.Ins,[src ':BURSt:NCYC 1']);
            fprintf(obj.Ins,[src ':BURSt:MODE TRIGgered']);
            fprintf(obj.Ins,['TRIGger' num2str(ch) ':SOURce ' TrigSource]); % why both???
            fprintf(obj.Ins, [':SOURce' num2str(ch) ':FUNCtion:PULSe:TRANsition:BOTH 0.000001']);
            OutputON(obj,ch);
            %         fprintf(obj.Ins,'OUTPut:TRIGger ON');
            %         fprintf(obj.Ins,['OUTPut:TRIGger:SOURce CH' num2str(ch)]);
            %         fprintf(obj.Ins,['TRIGger' num2str(ch) ':SOURce ' TrigSource]);
            %         fprintf(obj.Ins,['OUTPut:SYNC:SOURce CH' num2str(ch)]);
            %         if strcmpi(TrigSource,'BUS'); Trig(obj); end
            obj.readError(true);
        end
        
        function BurstSinV2(obj,ch,Ncyc,trig,freq,amp,offset,delay)
            if ~exist('ch','var') || isempty(ch); ch = 1; end
            if ~exist('offset','var') || isempty(offset); offset = 0; end
            if ~exist('trig','var') || isempty(trig); trig = 'BUS'; end
            if ~exist('delay','var') || isempty(delay); delay = 0; end
            src = [':SOURce' num2str(ch)];
            
            [~,wvfm] = obj.Output(ch,nan,true); if ~strcmp(wvfm{1}.name,'SIN'); fprintf(obj.Ins,[src ':FUNCtion SINusoid']); end
            
            fprintf(obj.Ins,[src sprintf(':FREQuency %g', freq)]);
            fprintf(obj.Ins,[src sprintf(':VOLTage %g', 2*amp)]);
            fprintf(obj.Ins,[src sprintf(':VOLTage:OFFSet %g', offset)]);
            
            fprintf(obj.Ins,[src ':BURSt:STATe ON']);
            fprintf(obj.Ins,['TRIGger' num2str(ch) ':DELay ' num2str(delay)]);
            fprintf(obj.Ins,[src ':BURSt:NCYC ' num2str(Ncyc)]);
            fprintf(obj.Ins,[src ':BURSt:MODE TRIGgered']);
            fprintf(obj.Ins,['TRIGger' num2str(ch) ':SOURce ' trig]);
            OutputON(obj,ch);
            %         fprintf(obj.Ins,'OUTPut:TRIGger ON');
            %         fprintf(obj.Ins,['OUTPut:TRIGger:SOURce CH' num2str(ch)]);
            %         fprintf(obj.Ins,['TRIGger' num2str(ch) ':SOURce ' TrigSource]);
            %         fprintf(obj.Ins,['OUTPut:SYNC:SOURce CH' num2str(ch)]);
            %         if strcmpi(TrigSource,'BUS'); Trig(obj); end;
            
            obj.readError(true);
        end
        
        function SingleTrapz(obj,ch,Vi,Vf,dt,TrigSource,dt0)
            if ~exist('ch','var') || isempty(ch) || isnan(ch); ch = 1; end
            if ~exist('TrigSource','var'); TrigSource = 'BUS'; end
            if ~exist('dt0','var') || isempty(dt0) || isnan(dt0); dt0 = 0; end
            src = [':SOURce' num2str(ch)];
            
            fprintf(obj.Ins,[src ':APPLy:RAMP ' num2str(1/(2*dt)) ',' num2str(abs(Vf-Vi)) ',' num2str((Vf+Vi)/2)]);
            fprintf(obj.Ins,[src ':BURSt:STATe ON']);
            fprintf(obj.Ins,['TRIGger' num2str(ch) ':SOURce BUS']);
            fprintf(obj.Ins,['TRIGger' num2str(ch) ':DELay ' num2str(dt0)]);
            fprintf(obj.Ins,[src ':BURSt:NCYC 1']);
            fprintf(obj.Ins,[src ':BURSt:MODE TRIGgered']);
            if Vi<=Vf; fprintf(obj.Ins,[src ':BURSt:PHASe +179.9']); end
            if Vf<Vi ; fprintf(obj.Ins,[src ':BURSt:PHASe -180']  ); end
            
            OutputON(obj,ch);
            
            fprintf(obj.Ins,'OUTPut:TRIGger ON');
            fprintf(obj.Ins,['OUTPut:TRIGger:SOURce CH' num2str(ch)]);
            fprintf(obj.Ins,['TRIGger' num2str(ch) ':SOURce ' TrigSource]);
            fprintf(obj.Ins,['OUTPut:SYNC:SOURce CH' num2str(ch)]);
            
            if strcmpi(TrigSource,'BUS'); Trig(obj); end
        end
        
        function LoadARB(obj,ch,Fs,data,Vpp,SigName)
            if ~exist('arbName','var'); SigName = obj.DefARBname; end
            if ~exist('ch','var') || isempty(ch) || isnan(ch); ch = 1; end
            if ~exist('Vpp','var') || isempty(Vpp) || isnan(Vpp); Vpp = 1; end
            %         buffer = length(data)*8;
            %         if get(obj.Ins,'OutputBufferSize')<(buffer+125)
            %             fclose(obj.Ins); set(obj.Ins,'OutputBufferSize',buffer+125); fopen(obj.Ins);
            %         end
            
            src = ['SOURce' num2str(ch)];
            data = (1*single(data(:)))/max(abs(data));
            
            fprintf(obj.Ins,[src ':DATA:VOLatile:CLEar']); %Clear volatile memory
            fprintf(obj.Ins,'FORM:BORD SWAP');  %configure the box to correctly accept the binary arb points
            Bytes = num2str(4*length(data)); % # of bytes
            header = [src ':DATA:ARBitrary ' SigName ', #' num2str(length(Bytes)) Bytes]; %create header
            binblockBytes = typecast(data.','uint8');  % convert datapoints to binary before sending
            fwrite(obj.Ins,[header binblockBytes], 'uint8'); % combine header and datapoints then send to instrument
            fprintf(obj.Ins,'*WAI'); % Make sure no other commands are exectued until arb is done downloading
            
            fprintf(obj.Ins,[src ':FUNCtion:ARBitrary ' SigName]); % set current arb waveform to defined arb vec
            fprintf(obj.Ins,['MMEM:STOR:DATA1 "INT:\' SigName '.arb"']); % store arb in intermal NV memory
            
            fprintf(obj.Ins,[src ':FUNCtion:ARB:SRATe ' num2str(Fs)]);%set sample rate
            fprintf(obj.Ins,[src ':FUNCtion ARB']); % turn on arb function
            fprintf(obj.Ins,[src ':VOLT ' num2str(Vpp)]); %send amplitude command
            %         fprintf(['Arb waveform downloaded to channel ' num2str(ch) '\n\n']); %print waveform has been downloaded
        end
        
        function ARB(obj,ch,delay,is_burst,offset,trigSource,Amp_pp)
            if ~exist('ch','var') || isempty(ch) || isnan(ch); ch = 1; end
            if ~exist('delay','var') || isempty(delay) || isnan(delay); delay = 0; end
            if ~exist('is_burst','var') || isempty(is_burst) || isnan(is_burst); is_burst = 1; end
            if ~exist('offset','var') || isempty(offset) || isnan(offset); offset = 0; end
            if ~exist('trigSource','var'); trigSource = 'EXTernal'; end
            
            obj.OutputOFF(ch);
            fprintf(obj.Ins,['TRIGger',num2str(ch),':SOURce ',trigSource]);
            fprintf(obj.Ins,['TRIGger',num2str(ch),':DELay ',num2str(delay)]);
            fprintf(obj.Ins,['SOUR',num2str(ch),':BURSt:STATe ',num2str(is_burst)']);
            fprintf(obj.Ins,['SOUR',num2str(ch),':VOLT:OFFS ',num2str(offset)]);
            if exist('Amp_pp','var'); fprintf(obj.Ins,['SOUR',num2str(ch),':VOLT ',num2str(Amp_pp)]); end
            obj.OutputON(ch);
        end
        
        % THIS FUNCITION IS OUTDATED! USE SinglePulse_v2 INSTEAD.
        function DelayedSinglePulse(obj,ch,V,dt,dt0)
            % THIS FUNCITION IS OUTDATED! USE SinglePulse_v2 INSTEAD.
            if ~exist('ch','var') || isempty(ch) || isnan(ch); ch = 1; end
            src = [':SOURce' num2str(ch)];
            
            fprintf(obj.Ins,[src ':APPLy:PULSe ' num2str(1/(2*(dt+dt0))) ',' num2str(2*V) ',' num2str(V)]);
            fprintf(obj.Ins,[src ':FUNCtion:PULSe:WIDTh ' num2str(dt)]);
            fprintf(obj.Ins,[src ':BURSt:STATe ON']);
            %         fprintf(obj.Ins,[src ':BURSt:PHASe ' num2str(360*(2*dt+dt0)/(2*(dt+dt0)))]);
            %         fprintf(obj.Ins,['TRIGger' num2str(ch) ':SOURce BUS']);
            fprintf(obj.Ins,['TRIGger' num2str(ch) ':SOURce EXTernal']);
            fprintf(obj.Ins,['TRIGger' num2str(ch) ':DELay ' num2str(dt0)]);
            fprintf(obj.Ins,[src ':BURSt:NCYC 1']);
            fprintf(obj.Ins,[src ':BURSt:MODE TRIGgered']);
            OutputON(obj,ch);
            
            %         fprintf(obj.Ins,'OUTPut:TRIGger ON');
            %         fprintf(obj.Ins,['OUTPut:TRIGger:SOURce CH' num2str(ch)]);
            %         fprintf(obj.Ins,['TRIGger' num2str(ch) ':SOURce BUS']);
            fprintf(obj.Ins,['OUTPut:SYNC:SOURce CH' num2str(ch)]);
            
            %         Trig(obj);
        end
        
        function state_out = Burst(obj,ch,state_in,source)
            if exist('state_in','var') && ~isempty(state_in)
                state_out = state_in;
                if state_in; BurstON(obj,ch); else BurstOFF(obj,ch); end
                if exist('source','var') && ischar(source); fprintf(obj.Ins,['TRIGger' num2str(ch) ':SOURce ' source]); end
            else state_out = num2str(query(obj.Ins, [':SOURce' num2str(ch) ':BURSt:STATe?']));
            end
        end
        
        function BurstOFF(obj,ch)
            if ~exist('ch','var') || isempty(ch); ch = 1; end
            fprintf(obj.Ins, [':SOURce' num2str(ch) ':BURSt:STATe OFF']);
        end
        
        function BurstON(obj,ch,source)
            if ~exist('ch','var') || isempty(ch); ch = 1; end
            fprintf(obj.Ins, [':SOURce' num2str(ch) ':BURSt:STATe ON']);
            if exist('source','var') && ischar(source); fprintf(obj.Ins,['TRIGger' num2str(ch) ':SOURce ' source]); end
        end
        
        function BurstON1(obj,ch,Nburst,source)
            if ~exist('ch','var') || isempty(ch); ch = 1; end
            if ~exist('Nburst','var') || isempty(Nburst) || ~isnumeric(Nburst) || isnan(Nburst); Nburst = 1; end
            fprintf(obj.Ins, [':SOURce' num2str(ch) ':BURSt:STATe ON']);
            %         fprintf(Obj.Ins,[':SOURce1:BURSt:NCYC INF']);
            fprintf(obj.Ins,[':SOURce' num2str(ch) ':BURSt:NCYC ' num2str(Nburst)]);
            if exist('source','var') && ischar(source); fprintf(obj.Ins,['TRIGger' num2str(ch) ':SOURce ' source]); end
        end
        
        function BurstON_INF(obj,ch)
            if ~exist('ch','var') || isempty(ch); ch = 1; end
            fprintf(obj.Ins, [':SOURce' num2str(ch) ':BURSt:STATe ON']);
            fprintf(Obj.Ins,':SOURce1:BURSt:NCYC INF');
        end
        
        function TrigOutOn(obj,ch)
            if ~exist('ch','var') || isempty(ch); ch = 1; end
            fprintf(obj.Ins, ':OUTPut:SYNC ON');
            fprintf(obj.Ins, [':OUTPut:SYNC:SOURce CH' num2str(ch) ]);
        end
        
        function polarity = Polarity(obj,ch_vec,sign)
            polarity = zeros(size(ch_vec));
            for ch = ch_vec
                temp = query(obj.Ins,['Output' num2str(ch) ':POL?']);
                polarity(ch==ch_vec) = strcmpi(temp(1:3),'NOR') - strcmpi(temp(1:3),'INV');
                if exist('sign','var') && ( (isnumeric(sign) && sign<0) || (ischar(sign) && strcmpi(sign(1:3),'INV')) )
                    fprintf(obj.Ins,['Output' num2str(ch) ':POL INV']);
                elseif exist('sign','var') && ( (isnumeric(sign) && sign>0) || (ischar(sign) && strcmpi(sign(1:4),'NORM')) )
                    fprintf(obj.Ins,['Output' num2str(ch) ':POL NORM']);
                end
            end
        end
        
        function VLim(obj,ch,vlim)
            if length(vlim)==1; vlim = vlim*[-1,1]; end
            fprintf(obj.Ins,['SOURce' num2str(ch) ':VOLTage:LIMit:LOW ' num2str(vlim(1))]);
            fprintf(obj.Ins,['SOURce' num2str(ch) ':VOLTage:LIMit:HIGH ' num2str(vlim(2))]);
            fprintf(obj.Ins,['SOURce' num2str(ch) ':VOLTage:LIMit:STATe ON']);
        end
        
        function delay = Delay(obj,ch,delay)
            if ~exist('ch','var'); ch = 1:2; end
            if ~exist('delay','var')
                delay = zeros(size(ch));
                for Ind=ch
                    delay(Ind==ch) = str2double(query(obj.Ins,['TRIGger' num2str(Ind) ':DELay?']));
                end
                return
            end
            if length(ch)~=length(delay); delay = delay*ones(size(ch)); end
            for Ind = ch
                fprintf(obj.Ins,['TRIGger' num2str(Ind) ':DELay ' num2str(delay(Ind==ch))]);
            end
        end
        
        function Abort(obj); fprintf(obj.Ins, 'Abort'); end
        
        function setTrig(obj,ch,mode); fprintf(obj.Ins, sprintf(':TRIGger%d:SOURce %s',ch,mode)); end % mode = IMM / EXT/ TIM / BUS
        
        function SetTrigger(obj,ch,mode); fprintf(obj.Ins, sprintf(':TRIGger%d:SOURce %s',ch,mode)); end % mode = IMM / EXT/ TIM / BUS
        
        function outtext = getChannelState(obj, ch)
            fprintf(obj.Ins,['SOURce' num2str(ch) ':VOLT?']);       vscale = fscanf(obj.Ins);
            fprintf(obj.Ins,['SOURce' num2str(ch) ':VOLT:OFFS?']);  voffset = fscanf(obj.Ins);
            fprintf(obj.Ins,['SOURce' num2str(ch) ':FUNC?']);       func = fscanf(obj.Ins);
            fprintf(obj.Ins,['SOURce' num2str(ch) ':FREQ?']);       freq = fscanf(obj.Ins);
            outstate = query(obj.Ins, ['OUTP' num2str(ch) '?']);
            
            outtext = ['CH: ' num2str(ch) ' ' 'Output' outstate ];
            outtext = [outtext ' ' 'Vscale: ' num2str(vscale)];
            outtext = [outtext ' ' 'Voffset: ' num2str(voffset)];
            outtext = [outtext ' ' 'Func: ' func];
            outtext = [outtext ' ' 'freq: ' freq];
        end
        
        function out = getChannelStateNum(obj, ch)
            % returns numeric values of instrument channel state
            fprintf(obj.Ins,['SOURce' num2str(ch) ':VOLT?']);       vscale = str2double(fscanf(obj.Ins));
            fprintf(obj.Ins,['SOURce' num2str(ch) ':VOLT:OFFS?']);  voffset = str2double(fscanf(obj.Ins));
            fprintf(obj.Ins,['SOURce' num2str(ch) ':FUNC?']);       func = fscanf(obj.Ins);
            fprintf(obj.Ins,['SOURce' num2str(ch) ':FREQ?']);       freq = str2double(fscanf(obj.Ins));
            outstate = str2double(query(obj.Ins, ['OUTP' num2str(ch) '?']));
            
            out.ch = ch;                % the channel in instrument
            out.state = outstate;       % binary state of channel: 1 = on, 0 = off
            out.vscale = vscale;        % voltage value in VOLTS
            out.voffset = voffset;      % voltage offset in VOLTS
            out.func = func(1:end - 1); % the waveform: 'SIN', 'PULSE', 'DC', ...
            out.freq = freq;            % the waveform frequency in Hz
        end
        
        function SweepConfig(obj,ch,type,timev,freq_startv,freq_stopv,hold_timev,return_timev)
            
            ch = num2str(ch);
            time = num2str(timev);
            freq_start = num2str(freq_startv);
            freq_stop = num2str(freq_stopv);
            hold_time = num2str(hold_timev);
            return_time=num2str(return_timev);
            
            % if ~exist('type','var') || isempty(type) || isnan(type); type='LIN'; end;
            if ~exist('type','var') || isempty(type) ; type='LIN'; end;
            if ~exist('time','var') || isempty(time) || isnan(timev); time = 1; end;
            if ~exist('hold_time','var') || isempty(hold_time) || isnan(hold_timev); hold_time=0; end;
            if ~exist('return_time','var') || isempty(return_time) || isnan(return_timev); return_time=0; end;
            
            fprintf(obj.Ins,['SOURce' ch ':SWEep:SPACing ' type]);
            fprintf(obj.Ins,['SOURce' ch ':SWEep:TIME ' time]);
            fprintf(obj.Ins,['SOURce' ch ':FREQuency:STARt ' freq_start]);
            fprintf(obj.Ins,['SOURce' ch ':FREQuency:STOP ' freq_stop]);
            fprintf(obj.Ins,['SOURce' ch ':SWEep:HTIMe ' hold_time]);
            fprintf(obj.Ins,['SOURce' ch ':SWEep:RTIMe ' return_time]);
        end
        
        function SweepON(obj,ch)
            ch = num2str(ch);
            fprintf(obj.Ins,['SOURce' ch ':SWEep:STATe ON']);
        end
        
        function SweepOFF(obj,ch)
            ch = num2str(ch);
            fprintf(obj.Ins,['SOURce' ch ':SWEep:STATe OFF']);
        end
        
        function setInetrnalFM(obj,ch,Amp,modFreq, modShape)
            if ~exist('ch','var') || isempty(ch) || isnan(ch); ch = 1; end
            if ~exist('modShape','var') || isempty(ch) || isnan(ch); modShape='SINusoid'; end
            fprintf(obj.Ins,[':SOURce' num2str(ch) ':FM:SOURce INTernal']);
            fprintf(obj.Ins,[':SOURce' num2str(ch) ':FM:INTernal:FUNCtion:SHAPe ' modShape]);
            fprintf(obj.Ins,[':SOURce' num2str(ch) ':FM:INTernal:FREQuency ' num2str(modFreq)]);
            fprintf(obj.Ins,[':SOURce' num2str(ch) ':FM:DEViation ' num2str(Amp)]);
        end
        
        function FMmodulationON(obj,ch)
            if ~exist('ch','var') || isempty(ch) || isnan(ch); ch = 1; end
            fprintf(obj.Ins,[':SOURce' num2str(ch) ':FM:STATe 1']);
        end
        
        function FMmodulationOFF(obj,ch)
            if ~exist('ch','var') || isempty(ch) || isnan(ch); ch = 1; end
            fprintf(obj.Ins,[':SOURce' num2str(ch) ':FM:STATe 0']);
        end
        
        function [stateOut, devOut] = FMmodulation(obj, ch, stateIn, modDev)
            % read / set FM modulation state and dev values on channel "ch"
            if ~exist('ch','var') || isempty(ch) || isnan(ch); ch = 1; end
            if exist('stateIn', 'var') && ~isempty(stateIn) && ~isnan(stateIn)
                fprintf(obj.Ins,[':SOURce' num2str(ch) ':FM:STATe '  num2str(stateIn)]);
                if exist('modDev', 'var') && ~isempty(modDev) && ~isnan(modDev)
                    fprintf(obj.Ins,[':SOURce' num2str(ch) ':FM:DEV '  num2str(modDev)]);
                end
            end
            if nargout >= 2
                stateOut = str2double(query(obj.Ins,['SOURce' num2str(ch) ':FM:STATe?']));
                devOut = str2double(query(obj.Ins,['SOURce' num2str(ch) ':FM:DEV?']));
            end
        end
        
        function inboundsflag = sigInbounds(obj,S,ch)
            R = obj.getVrange(ch);
            inboundsflag = ( max(S) < R*[0.005;0.995] ) && ( min(S) > R*[0.995;0.005] );
        end
        
        function [stateOut,isChSyncOut] = Syncio(obj, ch, stateIn)
            if ~exist('ch','var'); ch = 1; end
            ch_str = ['CH' num2str(ch)];
            if exist('stateIn', 'var') && ~isempty(stateIn) && ~isnan(stateIn)
                fprintf(obj.Ins, ['OUTPut:SYNC:SOURce ' ch_str]);
                fprintf(obj.Ins, [':OUTPut:SYNC:MODE NORMal']);
                fprintf(obj.Ins, [':OUTPut:SYNC:STATe ' num2str(stateIn)]);
            end
            if nargout >= 1
                stateOut = str2double(query(obj.Ins,[':OUTPut:SYNC:STATe?']));
                if stateOut 
                    chOut = query(obj.Ins,[':OUTPut:SYNC:SOURce?']);
                    isChSyncOut = double(str2double(chOut(3)) == ch);
                else
                    isChSyncOut = 0;
                end
            end
        end

    end
end

