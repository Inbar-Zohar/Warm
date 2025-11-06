% T = WW1071('192.168.1.26' , 'tcpip');
% 

classdef WW1071
% Updates:
% 16/05/19 - Roy: Added function ContinuousTrigOut (line 108)
% 04/11/19 - Roy: Fixed bug in SinglePulse in the definition of "dt" (line 85)
%            Roy: Fixed bug in Burst in RunMode: continuous vs triggered (line 90)
%            Roy: Added setSingleTrigOut function as a replacement for TrigOut (line 114)
%            Roy: Fixed bug in TrigOut in RunMode: continuous vs triggered (line 109)
% 11/02/20 - Roy: Added "fopen" in creator, in comment (line 28)
% 17/05/20 - Avraham Berrebi: add SinglePulseDelay (line 132)

   properties(Transient)
   Ins
   end
   methods
     function obj = WW1071(address,connection_type)
            if ~exist('address','var'); address = '192.168.1.26'; end
         if exist('connection_type','var') && strcmpi(connection_type,'GPIB')
           obj.Ins = gpib('NI', 0, address,'Timeout', 1.0);
        else
            obj.Ins = instrfind('Type', 'tcpip', 'RemoteHost', address, 'RemotePort', 23, 'Tag', '');
            if isempty(obj.Ins)
%                 obj.Ins = tcpip('192.168.1.30', 23, 'Timeout', 2);
                obj.Ins = tcpip(address, 23, 'Timeout', 2);
                set(obj.Ins, 'Name', ['TaborWW1071-' address]);
                else
                fclose(obj.Ins);
                obj.Ins = obj.Ins(1);
            end
         end
%          fopen(obj.Ins); % Gets open at the PID implementation
     end
     
      function data = IDN(obj); data = query(obj.Ins,'*IDN?'); end; %IDN
      function Reset(obj); fprintf(obj.Ins,'*RST'); end;
      function OutputON(obj);  fprintf(obj.Ins,'OUTP 1'); end;
      function OutputOFF(obj); fprintf(obj.Ins,'OUTP 0'); end;
      function [stt,wvfm] = Output(obj,sw,wvfm_flag)
          if exist('sw','var');
              stt = str2double(query(obj.Ins,'OUTPut?'));
              if exist('wvfm_flag','var') && wvfm_flag
                  temp = textscan(query(obj.Ins,'APPL?'),'%f,%f,%f');
                  wvfm.name = query(obj.Ins,':FUNCTion:SHAPe?');
                  while strcmp(wvfm.name(end),char(13)) || strcmp(wvfm.name(end),char(10))
                      wvfm.name(end) = [];
                  end
                  for J = 1:length(temp); wvfm.p(J)= temp{J}(1); end;
              end
          else fprintf(obj.Ins,sprintf('OUTPut %d',sw)); stt = sw; wvfm = [];
          end
      end
      
      function stt = OutputState(obj); stt = strcmpi(query(obj.Ins,'OUTP?'),sprintf('1\r\n')); end;
      function SyncON(obj);  fprintf(obj.Ins,'OUTPut:SYNC 1'); end;
      function SyncOFF(obj); fprintf(obj.Ins,'OUTPut:SYNC 0'); end;
      function stt = SyncOut(obj,sw)
          if ~exist('sw','var'); stt = str2double(query(obj.Ins,'OUTPut:SYNC?'));
          else fprintf(obj.Ins,sprintf('OUTPut:SYNC %d',sw)); stt = [];
          end
      end
      
      function Sin(obj,freq,amp,offs)
      data=[num2str(freq),',',num2str(amp),',',num2str(offs)];
      fprintf(obj.Ins,['APPL:SIN ',data]);
      end
      
      function Triangle(obj,freq,amp,offs)
      data=[num2str(freq),',',num2str(amp),',',num2str(offs)];
      fprintf(obj.Ins,['APPL:TRI ',data]);
      end
      
      function Ramp(obj,freq,amp,offs,delay,rise,fall)
      data=[num2str(freq),',',num2str(amp),',',num2str(offs),num2str(delay),',',num2str(rise),',',num2str(fall)];
      fprintf(obj.Ins,['APPL:RAMP ',data]);
      end
      
      function Pulse(obj,freq,amp,offs,delay,dt,rise,fall)
      data=[num2str(freq),',',num2str(amp),',',num2str(offs),',',num2str(delay),',',num2str(freq*dt*100),',',num2str(rise),',',num2str(fall)];
      fprintf(obj.Ins,['APPL:PULSe ',data]);
      end

      function Square(obj,freq,amp,offs,duty_cycle)
%       if ~exist('duty_cycle','var'); duty_cycle = 0.5; end;
      data=[num2str(freq),',',num2str(amp),',',num2str(0.5*offs),',',num2str(duty_cycle)];
      fprintf(obj.Ins,['APPL:SQU ',data]);
      end
      
      function DC(obj,amp)
      data = num2str(abs(amp));
      sgn = num2str(sign(amp)*100);
      fprintf(obj.Ins,['VOLT ',data]);
      fprintf(obj.Ins,['APPL:DC ',sgn]);
      end
      
      function SinglePulse(obj,dt,amp,trig_now)
      fprintf(obj.Ins,'INIT:CONT OFF');
      fprintf(obj.Ins,'TRIGger:BURSt ON');
      fprintf(obj.Ins,'TRIGger:COUNt 1');
      fprintf(obj.Ins,'TRIGger:SOURce:ADVance EXTernal');
      if ~exist('amp','var'); amp = 5; end;
      if ~exist('dt','var'); dt = 1e-6; end;
      freq = 1/(2*dt); offs = amp/2; delay = 0; rise = 0; fall = 0;
      obj.Pulse(freq,amp,offs,delay,0.5/freq,rise,fall);
      if exist('trig_now','var') && trig_now; obj.Trig; end;
      end
    
      function Burst(obj,sw,count)
      TrigTab.Continuous(0); % fprintf(obj.Ins,':TRIGger:MODE TRIGgered');
      fprintf(obj.Ins,sprintf(':TRIGger:BURSt %d',sw));
      fprintf(obj.Ins,':TRIGger:SOURce:ADVance EXTernal');
      if exist('count','var')
          fprintf(obj.Ins,[':TRIGger:BURSt:COUNt ' num2str(count)]);
      end
      end
      
      function Continuous(obj,sw)
      if ~exist('sw','var'); sw=1; end;
      fprintf(obj.Ins,['INIT:CONT ' num2str(sw)]);
      end
      
      function Trig(obj); fprintf(obj.Ins,'*TRG'); end;
      
      function TrigOut(obj,trig_now)
      if ~exist('trig_now','var'); trig_now = false; end;
      TrigTab.Continuous(0);
      obj.OutputON;
      obj.SinglePulse(10e-6,5,trig_now);
      end

      function setSingleTrigOut(obj,trig_now)
      if ~exist('trig_now','var'); trig_now = false; end;
      Continuous(obj,0);
      obj.OutputON;
      obj.SinglePulse(10e-6,5,trig_now);
      end

      function ContinuousTrigOut(obj,freq)
      if ~exist('freq','var'); freq = 1; end;
      Continuous(obj,1);
      Pulse(obj,freq,5,2.5,0,max(1e-3,0.2e-2/freq),0,0);
      OutputON(obj);
      end
      
      function SinglePulseDelay(obj,freq,width,amp,delay)
        fprintf(obj.Ins,sprintf('INST'));    
        fprintf(obj.Ins,'TRIGger:MODE:TRIG');
        fprintf(obj.Ins,'TRIGger:BURSt ON');
        fprintf(obj.Ins,'TRIGger:COUNt 1');
        fprintf(obj.Ins,'INIT:CONT ON');
        fprintf(obj.Ins,'INIT:CONT OFF');
        fprintf(obj.Ins,'TRIGger:SOURce:ADVance EXTernal');
        if ~exist('amp','var'); amp = 5; end;
        if ~exist('delay','var'); delay = 0; end;
        if width < 0.001 * 1/freq
            print('you need higher DDGfreq or freq, the minimum width is 0.1%')
        else
            width = width * freq * 100; % width = %up (1-99.99)
        end
        delay = delay * freq * 100;
        offs = amp/2; rise = 0; fall = 0; 
        data=[num2str(freq),',',num2str(abs(amp)),',',num2str(offs),',',num2str(delay),',',num2str(width),',',num2str(rise),',',num2str(fall)];

        fprintf(obj.Ins,['APPL:PULSe ',data]);
      end
      
   end
end