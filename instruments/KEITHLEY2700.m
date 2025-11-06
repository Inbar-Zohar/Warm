
% k = KEITHLEY2700(6,1);
% fopen(k.Ins)


% k.IDN
% k.set_auto_resolution;
% k.set_Filter('MOVING', 1, 10);
% k.set_rate(5);
% k.set_Resistance_meas;
% k.get_resistance


classdef KEITHLEY2700


   properties
   Ins
   end
   methods
       
      % Finding the address: 1) GPIB botton (+shift) and arrow up 
      % 2) At the PC -> With national instrument software
       function obj = KEITHLEY2700(address,is_visa)
           obj.Ins = gpib('NI', 0, address,'Timeout', 1.0);
           if exist('is_visa','var') && is_visa
               obj.Ins = visa('NI',['GPIB0::' num2str(address) '::INSTR']);
           end
       end
       
      function data = IDN(obj) 
          fprintf(obj.Ins,'*IDN?');
          data = fscanf(obj.Ins);
      end
      

      % KEITHLEY2700.set_Filter('MOVING', 1, 10)
      function set_Filter(obj,type, window_size, aquisition_num_in_Avg  ) % for ohm2 measurement specifically    
          fprintf(obj.Ins,['RESistance:AVERage:TCONtrol ' type]); % Set filter mode. Type = 'MOVING' or 'REPEAT'
          fprintf(obj.Ins,['RESistance:AVERage:WINDow ' num2str(window_size)]); % Set filter window in %. Use: 1
          fprintf(obj.Ins,['RESistance:AVERage:COUNt ' num2str(aquisition_num_in_Avg)]); % Set filter aquisition num. Use: 10
          fprintf(obj.Ins,'RESistance:AVERage:STATe 1'); % Activating the filter
      end

      
      function Filter_OFF(obj)
          fprintf(obj.Ins,'RESistance:AVERage:STATe 0'); % Activating the filter
      end
      
      
      function set_auto_resolution(obj)
          fprintf(obj.Ins,'RESistance:RANGe:AUTO 1'); % Setting aouto-range - gives the best resolution
      end
      
      
      function set_rate(obj, PLCs) %Set rate for ?2 in PLCs. Use: 5
         fprintf(obj.Ins,['RESistance:NPLCycles ' num2str(PLCs)]); 
      end
      

      function set_Resistance_meas(obj)
        fprintf(obj.Ins, 'INITiate:CONTinuous OFF');
        fprintf(obj.Ins, 'FUNC "RES"');
%         fprintf(obj.Ins, ':TRIGger:COUNt 1');
        fprintf(obj.Ins, ':READ?');
        fscanf(obj.Ins);
      end
      
      
      function data = get_resistance(obj)

         % Assuming that we activated: set_Resistance_meas
        fprintf(obj.Ins, ':READ?');
        
%          fprintf(obj.Ins,'MEASure:RESistance?');    % GIVES the right answer but take too long (2 seconds)
         tmp = fscanf(obj.Ins);
         % Use regular expression to extract resistance value
         tmp = regexp(tmp, '\+(\d+\.\d+E\+\d+)OHM', 'tokens');
         data = str2double(tmp{1}{1});
      end
      
      
        function get_error(obj) % prints the type of error (if exist)
         fprintf(obj.Ins,'STAT:QUE?');
         fscanf(obj.Ins)

      end
      
      
      
      
      
      
      
       
   end 
end