classdef SIM900
   properties
   Ins
   SAport
   PIDAport
   PIDBport
   end
    methods
      function obj = SIM900(com)
%         obj.Ins = instrfind('Type', 'serial', 'Port', 'COM8', 'Tag', '');
        obj.Ins = visa('AGILENT',['ASRL' num2str(com) '::INSTR']);
        obj.SAport = 1; % NMR summing amplifier port
        obj.PIDAport = 4; % Xe-129 (a) prot
        obj.PIDBport = 6; % Xe-131 (b) prot
        obj.Ins.Timeout = 2;
      end
      
      %port 4 - Xe131 
      %port 6 - Xe129 
      
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
      
      function data = ShowPropGain(obj,port) %Proportional Gain
          fprintf(obj.Ins,['conn ' num2str(port) ' ,''!' num2str(port) 'XYZ''']);
          fprintf(obj.Ins,'GAIN?');
          data=str2double(fscanf(obj.Ins));
          fprintf(obj.Ins,['!' num2str(port) 'XYZ']);
      end
      
      function data = SetPropGain(obj,port,gain) %set Proportional Gain
          fprintf(obj.Ins,['conn ' num2str(port) ' ,''!' num2str(port) 'XYZ''']);
          fprintf(obj.Ins,['GAIN ' num2str(gain)]);
          fprintf(obj.Ins,'GAIN?');
          data=str2double(fscanf(obj.Ins));
          fprintf(obj.Ins,['!' num2str(port) 'XYZ']);
      end
          
      function data = ShowIntGain(obj,port) %Integral Gain
          fprintf(obj.Ins,['conn ' num2str(port) ' ,''!' num2str(port) 'XYZ''']);
          fprintf(obj.Ins,'INTG?');
          data=str2double(fscanf(obj.Ins));
          fprintf(obj.Ins,['!' num2str(port) 'XYZ']);
      end
      
      function data = SetIntGain(obj,port,gain) %set Int Gain
          fprintf(obj.Ins,['conn ' num2str(port) ' ,''!' num2str(port) 'XYZ''']);
          fprintf(obj.Ins,['INTG' num2str(gain)]);
          fprintf(obj.Ins,'INTG?');
          data=str2double(fscanf(obj.Ins));
          fprintf(obj.Ins,['!' num2str(port) 'XYZ']);
      end
      
      function data = ShowDerGain(obj,port) %Derivative Gain
          fprintf(obj.Ins,['conn ' num2str(port) ' ,''!' num2str(port) 'XYZ''']);
          fprintf(obj.Ins,'DERV?');
          data=str2double(fscanf(obj.Ins));
          fprintf(obj.Ins,['!' num2str(port) 'XYZ']);
      end
      
      function data = SetDerGain(obj,port,gain) %set Int Gain
          fprintf(obj.Ins,['conn ' num2str(port) ' ,''!' num2str(port) 'XYZ''']);
          fprintf(obj.Ins,['DERV' num2str(gain)]);
          fprintf(obj.Ins,'DERV?');
          data=str2double(fscanf(obj.Ins));
          fprintf(obj.Ins,['!' num2str(port) 'XYZ']);
      end
      
      function data = SetDervGain(obj,port,gain) %set Derivativ Gain
          fprintf(obj.Ins,['conn ' num2str(port) ' ,''!' num2str(port) 'XYZ''']);
          fprintf(obj.Ins,['DERV' num2str(gain)]);
          fprintf(obj.Ins,'DERV?');
          data=str2double(fscanf(obj.Ins));
          fprintf(obj.Ins,['!' num2str(port) 'XYZ']);
      end
      
      function data = PropGain(obj,port,a) %Turn on/off Proportional Gain
          fprintf(obj.Ins,['conn ' num2str(port) ' ,''!' num2str(port) 'XYZ''']);
          if strcmp(a,'ON') || strcmp(a,'OFF')
            fprintf(obj.Ins,['PCTL ' a]);
          end
          fprintf(obj.Ins,'PCTL?');
          data=str2double(fscanf(obj.Ins));
          fprintf(obj.Ins,['!' num2str(port) 'XYZ']);
      end
      
      function data = IntGain(obj,port,a) %Turn on/off Integral Gain
          fprintf(obj.Ins,['conn ' num2str(port) ' ,''!' num2str(port) 'XYZ''']);
          if strcmp(a,'ON') || strcmp(a,'OFF')
            fprintf(obj.Ins,['ICTL ' a]);
          end
          fprintf(obj.Ins,'ICTL?');
          data=str2double(fscanf(obj.Ins));
          fprintf(obj.Ins,['!' num2str(port) 'XYZ']);
      end
      
      function data = DerGain(obj,port,a) %Turn on/off Derivative Gain
          fprintf(obj.Ins,['conn ' num2str(port) ' ,''!' num2str(port) 'XYZ''']);
          if strcmp(a,'ON') || strcmp(a,'OFF')
            fprintf(obj.Ins,['DCTL ' a]);
          end
          fprintf(obj.Ins,'DCTL?');
          data=str2double(fscanf(obj.Ins));
          fprintf(obj.Ins,['!' num2str(port) 'XYZ']);
      end
      
      
      %%%Controller configuration
      
      
        function data = UpLIM(obj,port,lim) %Controller Upper limit output [Volts]
          fprintf(obj.Ins,['conn ' num2str(port) ' ,''!' num2str(port) 'XYZ''']);
          fprintf(obj.Ins,['ULIM ' num2str(lim)]);
          fprintf(obj.Ins,'ULIM?');
          data=str2double(fscanf(obj.Ins));
          fprintf(obj.Ins,['!' num2str(port) 'XYZ']);
        end
      
      
        function data = LowLIM(obj,port,lim) %Controller Lower limit output [Volts]
          fprintf(obj.Ins,['conn ' num2str(port) ' ,''!' num2str(port) 'XYZ''']);
          fprintf(obj.Ins,['LLIM ' num2str(lim)]);
          fprintf(obj.Ins,'LLIM?');
          data=str2double(fscanf(obj.Ins));
          fprintf(obj.Ins,['!' num2str(port) 'XYZ']);
        end
      
      function data = getPIDA_IntGain(obj);  data = obj.ShowIntGain(obj.PIDAport);  end % returns Xe-129 Integral Gain
      function data = getPIDB_IntGain(obj);  data = obj.ShowIntGain(obj.PIDBport);  end % returns Xe-131 Integral Gain
      function data = getPIDA_PropGain(obj); data = obj.ShowPropGain(obj.PIDAport); end % returns Xe-129 Proportional Gain
      function data = getPIDB_PropGain(obj); data = obj.ShowPropGain(obj.PIDBport); end % returns Xe-131 Proportional Gain
      function data = getPIDA_DerGain(obj);  data = obj.ShowDerGain(obj.PIDAport);  end % returns Xe-129 Derivative Gain
      function data = getPIDB_DerGain(obj);  data = obj.ShowDerGain(obj.PIDBport);  end % returns Xe-131 Derivative Gain
      
      function data = setPIDA_IntGain(obj,g);  data = obj.SetIntGain(obj.PIDAport,g);  end % returns Xe-129 Integral Gain
      function data = setPIDB_IntGain(obj,g);  data = obj.SetIntGain(obj.PIDBport,g);  end % returns Xe-131 Integral Gain
      function data = setPIDA_PropGain(obj,g); data = obj.SetPropGain(obj.PIDAport,g); end % returns Xe-129 Proportional Gain
      function data = setPIDB_PropGain(obj,g); data = obj.SetPropGain(obj.PIDBport,g); end % returns Xe-131 Proportional Gain
      function data = setPIDA_DerGain(obj,g);  data = obj.SetDerGain(obj.PIDAport,g);  end % returns Xe-129 Derivative Gain
      function data = setPIDB_DerGain(obj,g);  data = obj.SetDerGain(obj.PIDBport,g);  end % returns Xe-131 Derivative Gain
      
      function data = Output(obj,port,state) %Integral Gain
          fprintf(obj.Ins,['conn ' num2str(port) ' ,''!' num2str(port) 'XYZ''']);
          if exist('state','var') && ~isempty(state) && ~isnan(state)
              if ischar(state)
                  switch lower(state); case 'm'; state = 0; case 'p'; state = 1; otherwise; error('ERROR! PID status is either 0 or 1'); end
              end
              state = double(state);
              fprintf(obj.Ins,['AMAN ' num2str(state)]);
          end
          data = str2double(query(obj.Ins,'AMAN?'));
          fprintf(obj.Ins,['!' num2str(port) 'XYZ']);
      end
      function data = PIDAOutput(obj,state)
          if ~exist('state','var'); state = nan; end
          data = obj.Output(obj.PIDAport,state);
      end
      function data = PIDBOutput(obj,state)
          if ~exist('state','var'); state = nan; end
          data = obj.Output(obj.PIDBport,state);
      end
    end
end


% PCTL(?) z 3 – 10 Proportional action ON/OFF
% ICTL(?) z 3 – 10 Integral action ON/OFF
% DCTL(?) z 3 – 10 Derivative action ON/OFF
% OCTL(?) z 3 – 10 Offset ON/OFF
% GAIN(?) {f} 3 – 10 Proportional Gain
% APOL(?) z 3 – 11 Controller Polarity
% INTG(?) {f} 3 – 11 Integral Gain
% DERV(?) {f} 3 – 11 Derivative Gain
% OFST(?) {f} 3 – 11 Output Offset
% Controller Configuration
% AMAN(?) z 3 – 12 Output (Manual Output/PID Control)
% INPT(?) z 3 – 12 Input (Internal/External Setpoint)
% SETP(?) {f} 3 – 12 New setpoint
% RAMP(?) z 3 – 12 Internal setpoint ramping ON/OFF
% RATE(?) {f} 3 – 12 Setpoint ramping Rate
% RMPS? 3 – 13 Setpoint ramping status
% STRT z 3 – 13 Pause or continue ramping
% MOUT(?) {f} 3 – 13 Manual Output
% ULIM(?) {f} 3 – 13 Upper Output Limit
% LLIM(?) {f} 3 – 14 Lower Output Limit