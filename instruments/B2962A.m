classdef B2962A
% Updates:
% 18/04/19 - Roy:
%	1- Changed defult trigger delay to 0 (line 196)
% 18/04/19 - Roy:
%	1- Added function OutputState (line 78)
%   2- Changed function Ramp with apply_immmedeate flag (line 162)
%   3- Added function Ramp3 (line 300)
% 01/03/20 - Roy:
%   1- heavily edited Pulse (line 90)
%   2- added output state to Output (line 84)
%   3- added function PulseUDEF (line 243)
%   4- added function readError (line 42)
%   5- added function Arm (line 216)
% 14/03/20 - Roy: added function VLim (line 354)
% 28/03/20 - Roy: 1- added function AdiabaticDCchange (?)
%                 2- added function ForceTrig (line 44)

   properties(Transient)
    Ins
   end
   methods
    function obj = B2962A(address)
        if ~exist('address','var'); address = '192.168.1.121'; end;
        obj.Ins = instrfind('Type','visa-tcpip','name',['TCPIP0::' address '::hislip0::INSTR']);
        if isempty(obj.Ins)
           obj.Ins = instrfind('Type','visa-tcpip','name',['TCPIP0::' address '::inst0::INSTR']);
        end
        if isempty(obj.Ins)
            obj.Ins = visa('AGILENT',['TCPIP0::' address '::hislip0::INSTR']);
        else
            fclose(obj.Ins);
            obj.Ins = obj.Ins(1);
        end
        set(obj.Ins,'InputBufferSize',2^10);
        set(obj.Ins,'OutputBufferSize',2^10);
        
        fopen(obj.Ins);
        fprintf(obj.Ins,':SOURce1:FUNCtion:MODE CURRent'); %% sets current mode
        fprintf(obj.Ins,':SOURce2:FUNCtion:MODE CURRent'); %% sets current mode
        fprintf(obj.Ins,':SOURce1:CURRent:RANGe:AUTO 1'); %% sets the range automatically
        fprintf(obj.Ins,':SOURce2:CURRent:RANGe:AUTO 1'); %% sets the range automatically
        fprintf(obj.Ins,':SENSe1:VOLTage:DC:PROTection:LEVel:BOTH 2'); %% limits the maximal voltage of the supplier (2V ~ 1A)
        fprintf(obj.Ins,':SENSe2:VOLTage:DC:PROTection:LEVel:BOTH 2'); %% limits the maximal voltage of the supplier (2V ~ 1A)
%         fclose(obj.Ins);
    end
    function ForceTrig(obj); fprintf(obj.Ins, '*TRG'); end;
    
    function data = IDN(obj); data = query(obj.Ins,'*IDN?'); end %IDN
    function [errCode,errString] = readError(obj,disp_flag)
        err = textscan(query(obj.Ins, ':SYSTem:ERRor?'),'%d%q','Delimiter',',');
        errCode = err{1}; errString = err{2}{1};
        if strcmp(errString(1:8),'No error')==1; errCode = false; errString = [];
        elseif exist('disp_flag','var') && disp_flag; disp(errString);
        end
    end
    
    function Current(obj,src,I)
        if ~exist('src','var'); src = 1; end;
        fprintf(obj.Ins,[':SOURce' num2str(src) ':CURRent:LEVel:IMMediate:AMPLitude ' num2str(I)]);
    end
    function data = ReadCurrent(obj,src)
        if ~exist('src','var'); src = 1; end;
        data = str2double(query(obj.Ins,...
            [':SOURce' num2str(src) ':CURRent:LEVel:IMMediate:AMPLitude?']));
    end
        
    function OutputON(obj,src)
        if ~exist('src','var'); src = 1; end;
        fprintf(obj.Ins,[':OUTPut' num2str(src) ':STATe 1']);
    end
    function OutputOFF(obj,src)
        if ~exist('src','var'); src = 1; end;
        fprintf(obj.Ins,[':OUTPut' num2str(src) ':STATe 0']);
    end
    function stt = State(obj,src)
        if ~exist('src','var'); src = 1:2; end;
        stt = zeros(size(src));
        for s = src
            stt(s) = str2double(query(obj.Ins,...
                [':OUTPut' num2str(s) ':STATe?']));
        end
        if length(src)==1 && src==2; stt = stt(2); end;
    end
	function stt = OutputState(obj,src)
        if ~exist('src','var'); src = 1:2; end;
        stt = zeros(size(src));
        for s = src
            stt(s) = str2double(query(obj.Ins,...
                [':OUTPut' num2str(s) ':STATe?']));
        end
        if length(src)==1 && src==2; stt = stt(2); end;
    end
    function stt = Output(obj,ch,stt)
        if ~exist('ch','var'); ch = 1:2; end;
        if ~exist('stt','var')
            for I=ch; stt(I==ch) = OutputState(obj,ch); end;
        else
            for I=ch; if stt(I==ch); OutputON(I); else OutputOFF(I); end; end;
        end
    end
    
    function Pulse(obj,src,I,dt,freq,I0,Npuls,dt0,isArmFlag)
        % src   = ch# to apply the output for (1 or 2?)
        % I     = current (in Amp) during the pulse
        % dt    = pulse duration
        % freq  = pulse repitition rate (PRR)
        % I0    = current (in Amp) to apply after the pulse
        % Npuls = number of pulses to output (they are outputed in PRR?)
        % dt0   = time delay between trigger and pulse output
        if ~exist('Npuls','var') || isempty(Npuls) || (isnumeric(Npuls) && isnan(Npuls)); Npuls=1; end;
%         if ~exist('Npuls','var') || isempty(Npuls) || isnan(Npuls); Npuls='INFinity'; end;
        if ~exist('I0','var') || isempty(I0) || isnan(I0); I0 = 0; end;
        if ~exist('dt0','var') || isempty(dt0) || isnan(dt0); dt0 = 0; end;
        if ~exist('freq','var') || isempty(freq) || isnan(freq); freq = 1/(2*dt); end;
        
        % cancel armmed sequenece (it causes problems if it exists)
        obj.Abort;
        
        % define DC value:
        fprintf(obj.Ins,sprintf(':SOURce%d:CURRent:LEVel:IMMediate:AMPLitude %g',src,I0)); %% changes immediately the current state to the value in Ampere
%         fprintf(obj.Ins,sprintf(':OUTPut%d:STATe %d',ch, 1));%% Turn on channel ch
        
        % define pulse:
        fprintf(obj.Ins,sprintf(':SOURce%d:CURRent:MODE %s',src,'ARB'));
        fprintf(obj.Ins,sprintf(':SOURce%d:ARB:FUNCtion:SHAPe %s',src,'SQUare'));
        if isnumeric(Npuls); Npuls = upper(num2str(Npuls)); end;
            fprintf(obj.Ins,sprintf(':SOURce%d:ARB:COUNt %s',src, Npuls));
        fprintf(obj.Ins,sprintf(':SOURce%d:ARB:CURRent:SQUare:STARt:LEVel %g',src,I0));
        fprintf(obj.Ins,sprintf(':SOURce%d:ARB:CURRent:SQUare:STARt:TIME %g',src,dt0));
        fprintf(obj.Ins,sprintf(':SOURce%d:ARB:CURRent:SQUare:TOP:LEVel %g',src,I));
        fprintf(obj.Ins,sprintf(':SOURce%d:ARB:CURRent:SQUare:TOP:TIME %g',src,dt));
        fprintf(obj.Ins,sprintf(':SOURce%d:ARB:CURRent:SQUare:END:TIME %g',src,1/freq-dt));
        
        if exist('isArmFlag','var') && isArmFlag; obj.Arm(src,I0); end
    end
    
    function Sin(obj,src,Amp,freq,I0,DC,Npuls)
        if ~exist('Npuls','var') || isempty(Npuls) || isnan(Npuls); Npuls='INFinity'; end; % Npuls=1;
        if ~exist('DC','var') || isempty(DC) || isnan(DC); DC = I0; end;
        % define DC value:
        if ~exist('I0','var') || isempty(I0) || isnan(I0); I0 = 0; end;
        fprintf(obj.Ins,sprintf(':SOURce%d:CURRent:LEVel:IMMediate:AMPLitude %g',src,I0)); %% changes immediately the current state to the value in Ampere
%         fprintf(obj.Ins,sprintf(':OUTPut%d:STATe %d',ch,1));%% Turn on channel ch
        
        % define sine:
        fprintf(obj.Ins,sprintf(':SOURce%d:CURRent:MODE %s',src,'ARB'));
        fprintf(obj.Ins,sprintf(':SOURce%d:ARB:FUNCtion:SHAPe %s',src,'SINusoid'));
        if isnumeric(Npuls); Npuls = num2str(Npuls); end;
            fprintf(obj.Ins,sprintf(':SOURce%d:ARB:COUNt %s',src, Npuls));
        fprintf(obj.Ins,sprintf(':SOURce%d:ARB:CURRent:SINusoid:AMPLitude %g',src,Amp));
        fprintf(obj.Ins,sprintf(':SOURce%d:ARB:CURRent:SINusoid:FREQuency %g',src,freq));
        fprintf(obj.Ins,sprintf(':SOURce%d:ARB:CURRent:SINusoid:OFFSet %g',src,DC));
    end
    
    function Ramp(obj,src,t0,I0,Rt,Amp,TOPt,Npuls,appl_immediate)
        % defined in B2961-90030, page 2-4 (55)
        if ~exist('Npuls','var') || isempty(Npuls) || isnan(Npuls); Npuls='INFinity'; end; % Npuls=1;
        if ~exist( 'I0' ,'var') || isempty( I0 ) || ~isnumeric( I0 ) || isnan( I0 );  I0  = 0; end;
        if ~exist( 't0' ,'var') || isempty( t0 ) || ~isnumeric( t0 ) || isnan( t0 );  t0  = 0; end;
        
        if ~exist('appl_immediate','var') || appl_immediate
            fprintf(obj.Ins,sprintf(':SOURce%d:CURRent:LEVel:IMMediate:AMPLitude %g',src,I0)); % changes immediately the current state to the value in Ampere
        end
%         fprintf(obj.Ins,sprintf(':OUTPut%d:STATe %d',ch,1));%% Turn on channel ch
        
        % define ramp: (page 22 in manual: literature.cdn.keysight.com/litweb/pdf/B2961-90020.pdf)
        % define ramp: in B2961-90030.pdf, page 2-4 (55)
        fprintf(obj.Ins,sprintf(':SOURce%d:CURRent:MODE %s',src,'ARB'));
        fprintf(obj.Ins,sprintf(':SOURce%d:ARB:FUNCtion:SHAPe %s',src,'RAMP'));
        if isnumeric(Npuls); Npuls = num2str(Npuls); end;
            fprintf(obj.Ins,sprintf(':SOURce%d:ARB:COUNt %s',src, Npuls));
        fprintf(obj.Ins,sprintf(':SOURce%d:ARB:CURRent:RAMP:STARt:TIME %g',src,t0));
        fprintf(obj.Ins,sprintf(':SOURce%d:ARB:CURRent:RAMP:STARt %g',src,I0));
        fprintf(obj.Ins,sprintf(':SOURce%d:ARB:CURRent:RAMP:RTIM %g',src,Rt));
        fprintf(obj.Ins,sprintf(':SOURce%d:ARB:CURRent:RAMP:END %g',src,I0+Amp));
        fprintf(obj.Ins,sprintf(':SOURce%d:ARB:CURRent:RAMP:END:TIME %g',src,TOPt));
    end
    function Triangle(obj,src,t0,I0,Rt,Amp,Ft,tend,Npuls)
        if ~exist('Npuls','var') || isempty(Npuls) || isnan(Npuls); Npuls='INFinity'; end; % Npuls=1;
        if ~exist( 'I0' ,'var') || isempty( I0 ) || ~isnumeric( I0 ) || isnan( I0 );  I0  = 0; end;
        if ~exist( 't0' ,'var') || isempty( t0 ) || ~isnumeric( t0 ) || isnan( t0 );  t0  = 0; end;
        if ~exist('tend','var') || isempty(tend) || ~isnumeric(tend) || isnan(tend); tend = 0; end;
        fprintf(obj.Ins,sprintf(':SOURce%d:CURRent:LEVel:IMMediate:AMPLitude %g',src,I0)); %% changes immediately the current state to the value in Ampere
%         fprintf(obj.Ins,sprintf(':OUTPut%d:STATe %d',ch,1));%% Turn on channel ch
        
        % define triangle: (page 22 in manual: literature.cdn.keysight.com/litweb/pdf/B2961-90020.pdf
        fprintf(obj.Ins,sprintf(':SOURce%d:CURRent:MODE %s',src,'ARB'));
        fprintf(obj.Ins,sprintf(':SOURce%d:ARB:FUNCtion:SHAPe %s',src,'TRIangle'));
        if isnumeric(Npuls); Npuls = num2str(Npuls); end;
            fprintf(obj.Ins,sprintf(':SOURce%d:ARB:COUNt %s',src, Npuls));
        fprintf(obj.Ins,sprintf(':SOURce%d:ARB:CURRent:TRIangle:STARt:TIME %g',src,t0));
        fprintf(obj.Ins,sprintf(':SOURce%d:ARB:CURRent:TRIangle:STARt %g',src,I0));
        fprintf(obj.Ins,sprintf(':SOURce%d:ARB:CURRent:TRIangle:RTIM %g',src,Rt));
        fprintf(obj.Ins,sprintf(':SOURce%d:ARB:CURRent:TRIangle:TOP %g',src,I0+Amp));
        fprintf(obj.Ins,sprintf(':SOURce%d:ARB:CURRent:TRIangle:FTIM %g',src,Ft));
        fprintf(obj.Ins,sprintf(':SOURce%d:ARB:CURRent:TRIangle:END:TIME %g',src,tend));
    end
    
    function Trapz(obj,src,t0,I0,Rt,Amp,TOPt,Ft,tend,Npuls)
        if ~exist('Npuls','var') || isempty(Npuls) || isnan(Npuls) || isinf(Npuls); Npuls='INFinity'; end; % Npuls=1;
        if ~exist( 'I0' ,'var') || isempty( I0 ) || ~isnumeric( I0 ) || isnan( I0 );  I0  = 0; end;
        if ~exist( 't0' ,'var') || isempty( t0 ) || ~isnumeric( t0 ) || isnan( t0 );  t0  = 0; end;
        if ~exist('tend','var') || isempty(tend) || ~isnumeric(tend) || isnan(tend); tend = 0; end;
        fprintf(obj.Ins,sprintf(':SOURce%d:CURRent:LEVel:IMMediate:AMPLitude %g',src,I0)); %% changes immediately the current state to the value in Ampere
%         fprintf(obj.Ins,sprintf(':OUTPut%d:STATe %d',ch,1));%% Turn on channel ch
        
        % define triangle: (page 22 in manual: literature.cdn.keysight.com/litweb/pdf/B2961-90020.pdf
        fprintf(obj.Ins,sprintf(':SOURce%d:CURRent:MODE %s',src,'ARB'));
        fprintf(obj.Ins,sprintf(':SOURce%d:ARB:FUNCtion:SHAPe %s',src,'TRAP'));
        if isnumeric(Npuls); Npuls = num2str(Npuls); end;
            fprintf(obj.Ins,sprintf(':SOURce%d:ARB:COUNt %s',src, Npuls));
        fprintf(obj.Ins,sprintf(':SOURce%d:ARB:CURRent:TRAP:STARt:TIME %g',src,t0));
        fprintf(obj.Ins,sprintf(':SOURce%d:ARB:CURRent:TRAP:STARt %g',src,I0));
        fprintf(obj.Ins,sprintf(':SOURce%d:ARB:CURRent:TRAP:RTIM %g',src,Rt));
        fprintf(obj.Ins,sprintf(':SOURce%d:ARB:CURRent:TRAP:TOP %g',src,Amp));
        fprintf(obj.Ins,sprintf(':SOURce%d:ARB:CURRent:TRAP:TOP:TIME %g',src,TOPt));
        fprintf(obj.Ins,sprintf(':SOURce%d:ARB:CURRent:TRAP:FTIM %g',src,Ft));
        fprintf(obj.Ins,sprintf(':SOURce%d:ARB:CURRent:TRAP:END:TIME %g',src,tend));
    end
    
    function Arm(obj,src,I0,delay)
        if ~exist('src','var') || isempty(src) || isnan(src); src = 1; end;
        % until trigger is on, output is in a DC mode on the offset value
        %     once triggr arrives it initializes the armmed waveform
        obj.Abort;
        if exist('I0','var') && ~isempty(I0) && ~isnan(I0) && isnumeric(I0)
            fprintf(obj.Ins,sprintf(':SOURce%d:CURRent:LEVel:IMMediate:AMPLitude %g',src,I0)); % changes immediately the current state to I0 value (in Ampere)
        end
        if ~obj.Output(src); fprintf(obj.Ins,sprintf(':OUTPut%d:STATe %d',src,1)); end; %% Turn on channel ch
        fprintf(obj.Ins,sprintf(':TRIGger%d:TRANsient:SOURce:SIGNal %s',src, 'EXT9'));
        fprintf(obj.Ins,':ARM(@%d):TRANsient:LAYer:IMMediate',src);
        if ~exist('delay','var'); delay = 0; end;
        fprintf(obj.Ins,[':TRIG:TRANsient:DELay ' num2str(delay)]);
        fprintf(obj.Ins,':INITiate:IMMediate:ALL');
    end
    
    function Trig(obj,src,I0,delay)
        if ~exist('src','var') || isempty(src) || isnan(src); src = 1; end;
        % until trigger is on, output is in a DC mode on the offset value
        %     once triggr arrives it initializes the armmed waveform
        if exist('I0','var') && ~isempty(I0) && ~isnan(I0) && isnumeric(I0)
            fprintf(obj.Ins,sprintf(':SOURce%d:CURRent:LEVel:IMMediate:AMPLitude %g',src,I0)); % changes immediately the current state to I0 value (in Ampere)
        end
        fprintf(obj.Ins,sprintf(':OUTPut%d:STATe %d',src,1));%% Turn on channel ch
        fprintf(obj.Ins,sprintf(':TRIGger%d:TRANsient:SOURce:SIGNal %s',src, 'EXT9'));
        fprintf(obj.Ins,':ARM(@%d):TRANsient:LAYer:IMMediate',src);
        if ~exist('delay','var'); delay = 0; end;
        fprintf(obj.Ins,[':TRIG:ALL:DELay ' num2str(delay)]);
        fprintf(obj.Ins,':INITiate:IMMediate:ALL');
    end
    
    function Abort(obj); fprintf(obj.Ins,':ABORt:TRANsient (@1,2)'); end;
    
    function PulseUDEF(obj,src,I_vec,dt,Npuls,I0)
        if ~exist('Npuls','var') || isempty(Npuls) || (isnumeric(Npuls) && isnan(Npuls)); Npuls=1; end;
        
        obj.Abort;
        
        % define DC value:
        if ~exist('I0','var') || isempty(I0) || isnan(I0); I0 = I_vec(1); end;
        fprintf(obj.Ins,sprintf(':SOURce%d:CURRent:LEVel:IMMediate:AMPLitude %g',src,I0)); %% changes immediately the current state to the value in Ampere
        
        % define pulse:
        fprintf(obj.Ins,sprintf(':SOURce%d:CURRent:MODE %s',src,'ARB'));
        fprintf(obj.Ins,sprintf(':SOURce%d:ARB:FUNCtion:SHAPe %s',src,'UDEF'));
        if isnumeric(Npuls); Npuls = num2str(Npuls); end;
            fprintf(obj.Ins,sprintf(':SOURce%d:ARB:COUNt %s',src, Npuls));
        
        I_str = [];
        for Ind = 1:length(I_vec)
            I_str = [I_str num2str(I_vec(Ind))]; %#ok<AGROW>
            if Ind~=length(I_vec); I_str = [I_str ',']; end %#ok<AGROW>
        end
        fprintf(obj.Ins,sprintf(':SOURce%d:ARB:CURRent:UDEfined:TIME %g',src,dt));
        fprintf(obj.Ins,sprintf(':SOURce%d:ARB:CURRent:UDEfined:LEVel %s',src,I_str));
        
        if exist('isArmFlag','var') && isArmFlag; obj.Arm(src,I0); end
    end
    
    function Pulse2(obj,src,t1,I1,t2,I2,Npuls)
        if ~exist('Npuls','var') || isempty(Npuls) || (isnumeric(Npuls) && isnan(Npuls)); Npuls=1; end;
        
        % define DC value:
        fprintf(obj.Ins,sprintf(':SOURce%d:CURRent:LEVel:IMMediate:AMPLitude %g',src,I1)); %% changes immediately the current state to the value in Ampere
%         if strcmpi(query(obj.Ins,sprintf(':OUTPut%d:STATe?',src)),'off')
%             fprintf(obj.Ins,sprintf(':OUTPut%d:STATe %d',src,1));
%         end
        
        % define pulse:
        fprintf(obj.Ins,sprintf(':SOURce%d:CURRent:MODE %s',src,'ARB'));
        fprintf(obj.Ins,sprintf(':SOURce%d:ARB:FUNCtion:SHAPe %s',src,'SQUare'));
        if isnumeric(Npuls); Npuls = num2str(Npuls); end;
            fprintf(obj.Ins,sprintf(':SOURce%d:ARB:COUNt %s',src, Npuls));
            
        fprintf(obj.Ins,sprintf(':SOURce%d:ARB:CURRent:SQUare:STARt:TIME %g',src,t1));
        fprintf(obj.Ins,sprintf(':SOURce%d:ARB:CURRent:SQUare:STARt:LEVel %g',src,I1));
        fprintf(obj.Ins,sprintf(':SOURce%d:ARB:CURRent:SQUare:TOP:TIME %g',src,t2));
        fprintf(obj.Ins,sprintf(':SOURce%d:ARB:CURRent:SQUare:TOP:LEVel %g',src,I2));
        fprintf(obj.Ins,sprintf(':SOURce%d:ARB:CURRent:SQUare:END:TIME %g',src,0));
    end
    function Pulse3(obj,src,I1,t2,I2,t3,I3,t4_for_I2,Npuls)
        if ~exist('Npuls','var') || isempty(Npuls) || (isnumeric(Npuls) && isnan(Npuls)); Npuls=1; end;
        
        % define DC value:
        fprintf(obj.Ins,sprintf(':SOURce%d:CURRent:LEVel:IMMediate:AMPLitude %g',src,I1)); %% changes immediately the current state to the value in Ampere

        
        % define pulse:
        fprintf(obj.Ins,sprintf(':SOURce%d:CURRent:MODE %s',src,'ARB'));
        fprintf(obj.Ins,sprintf(':SOURce%d:ARB:FUNCtion:SHAPe %s',src,'SQUare'));
        if isnumeric(Npuls); Npuls = num2str(Npuls); end;
            fprintf(obj.Ins,sprintf(':SOURce%d:ARB:COUNt %s',src, Npuls));
            
        fprintf(obj.Ins,sprintf(':SOURce%d:ARB:CURRent:SQUare:STARt:TIME %g',src,t2));
        fprintf(obj.Ins,sprintf(':SOURce%d:ARB:CURRent:SQUare:STARt:LEVel %g',src,I2));
        fprintf(obj.Ins,sprintf(':SOURce%d:ARB:CURRent:SQUare:TOP:TIME %g',src,t3));
        fprintf(obj.Ins,sprintf(':SOURce%d:ARB:CURRent:SQUare:TOP:LEVel %g',src,I3));
        if ~exist('t4_for_I2','var') || isempty(t4_for_I2) || ~isnumeric(t4_for_I2) || isnan(t4_for_I2) || isinf(t4_for_I2)
            t4_for_I2 = 0;
        end
        fprintf(obj.Ins,sprintf(':SOURce%d:ARB:CURRent:SQUare:END:TIME %g',src,t4_for_I2));
        
        if strcmpi(query(obj.Ins,sprintf(':OUTPut%d:STATe?',src)),'off')
            fprintf(obj.Ins,sprintf(':OUTPut%d:STATe %d',src,1));
        end
        obj.Trig(src,I1);
    end
	
    function Ramp3(obj,src,I0,t1,I1,Rt,I2,TOPt,Npuls)
        % defined in B2961-90030, page 2-4 (55)
        if ~exist('Npuls','var') || isempty(Npuls) || isnan(Npuls); Npuls='INFinity'; end; % Npuls=1;
        if ~exist( 'I0' ,'var') || isempty( I0 ) || ~isnumeric( I0 ) || isnan( I0 );  I0  = 0; end;
        if ~exist( 'I1' ,'var') || isempty( I1 ) || ~isnumeric( I1 ) || isnan( I1 );  I1  = 0; end;
        if ~exist( 't1' ,'var') || isempty( t1 ) || ~isnumeric( t1 ) || isnan( t1 );  t1  = 0; end;
        
        fprintf(obj.Ins,sprintf(':SOURce%d:CURRent:LEVel:IMMediate:AMPLitude %g',src,I0)); % changes immediately the current state to the value in Ampere
%         fprintf(obj.Ins,sprintf(':OUTPut%d:STATe %d',ch,1));%% Turn on channel ch
        
        % define ramp: (page 22 in manual: literature.cdn.keysight.com/litweb/pdf/B2961-90020.pdf)
        % define ramp: in B2961-90030.pdf, page 2-4 (55)
        fprintf(obj.Ins,sprintf(':SOURce%d:CURRent:MODE %s',src,'ARB'));
        fprintf(obj.Ins,sprintf(':SOURce%d:ARB:FUNCtion:SHAPe %s',src,'RAMP'));
        if isnumeric(Npuls); Npuls = num2str(Npuls); end;
            fprintf(obj.Ins,sprintf(':SOURce%d:ARB:COUNt %s',src, Npuls));
        fprintf(obj.Ins,sprintf(':SOURce%d:ARB:CURRent:RAMP:STARt:TIME %g',src,t1));
        fprintf(obj.Ins,sprintf(':SOURce%d:ARB:CURRent:RAMP:STARt %g',src,I1));
        fprintf(obj.Ins,sprintf(':SOURce%d:ARB:CURRent:RAMP:RTIM %g',src,Rt));
        fprintf(obj.Ins,sprintf(':SOURce%d:ARB:CURRent:RAMP:END %g',src,I2));
        fprintf(obj.Ins,sprintf(':SOURce%d:ARB:CURRent:RAMP:END:TIME %g',src,TOPt));
    end
    
    function getLim = VLim(obj,setLim,ch)
        if ~exist('ch','var') || isempty(ch) || any(isnan(ch)); ch = 1:2; end
        if ~exist('setLim','var') || isempty(setLim) || any(isnan(setLim))
            for c=ch
                getLim(c==ch) = str2double(query(obj.Ins,[':SENSe' num2str(c) ...
                                ':VOLTage:DC:PROTection:LEVel?']));
            end
            return
        else
            if length(setLim)<length(ch); setLim = setLim*ones(size(ch)); end
            for c=ch
                fprintf(obj.Ins,[':SENSe' num2str(c) ':VOLTage:DC:'...
                    'PROTection:LEVel:BOTH ' num2str(setLim(c==ch))]);
            end
            getLim = setLim;
        end
    end
    
    function getLim = ILim(obj,setLim,ch)
        if ~exist('ch','var') || isempty(ch) || any(isnan(ch)); ch = 1:2; end
        if ~exist('setLim','var') || isempty(setLim) || any(isnan(setLim))
            for c=ch
                getLim(c==ch) = str2double(query(obj.Ins,[':SENSe' num2str(c) ...
                                ':CURRent:DC:PROTection:LEVel?']));
            end
            return
        else
            if length(setLim)<length(ch); setLim = setLim*ones(size(ch)); end
            for c=ch
                fprintf(obj.Ins,[':SENSe' num2str(c) ':CURRent:DC:'...
                    'PROTection:LEVel:BOTH ' num2str(setLim(c==ch))]);
            end
            getLim = setLim;
        end
    end
    
%     function AdiabaticDCchange(obj,src,Iend,Istart,dt,applyNowFlag,wT_instead_of_dt_flag)
%         if ~exist('Istart','var') || isempty(Istart) || ~isnumeric(Istart) || isnan(Istart)
%             obj.Abort; Istart = ReadCurrent(obj,src);
%         end
%         if ~exist('dt','var') || isempty(dt) || ~isnumeric(dt) || isnan(dt)
%             dw_dI = 1e7;
%             if abs(dw_dI*Istart)>1; dt = 2*pi/abs(dw_dI*Istart)*
   end
end