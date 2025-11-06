classdef BK9174
   properties(Transient)
   Ins
   end
   properties
   storage
   end
   methods
     function obj = BK9174(com)
     obj.Ins = serial(sprintf('COM%d',com), 'BaudRate', 57600, 'Terminator', {'CR/LF','CR/LF'});%Set serial port
     end
     
     function data = IDN(obj) %IDN
      fprintf(obj.Ins,'*IDN?');
      data = fscanf(obj.Ins);
     end
        
      function Reset(obj)
      fprintf(obj.Ins,'*RST');
      end
      function SetVolt(obj,ch,volt)
      if ch==1; fprintf(obj.Ins,sprintf(['SOUR:VOLT' ' %f'],volt)); 
      else fprintf(obj.Ins,sprintf(['SOUR:VOLT2' ' %f'],volt)); end;
      end
      
      function data = Volt(obj,ch)
      if ch==1; fprintf(obj.Ins,'VOLT?');
      elseif ch==2; fprintf(obj.Ins,'VOLT2?');
      end
      data=str2double(fscanf(obj.Ins));
      end
      
      function SetCurr(obj,ch,curr)
      if ch==1; fprintf(obj.Ins,sprintf(['SOUR:CURR' ' %f'],curr)); 
      elseif ch==2; fprintf(obj.Ins,sprintf(['SOUR:CURR2' ' %f'],curr));
      end
      end
      
      function data = Curr(obj,ch)
      if ch==1; fprintf(obj.Ins,'CURR?');
      elseif ch==2; fprintf(obj.Ins,'CURR2?');
      end
      data=str2double(fscanf(obj.Ins));
      end
      
      function Output(obj,ch,sw)
      fprintf(obj.Ins,sprintf(['OUT' num2str(ch) ' %f'],sw));
      end
      
      function LimitVolt(obj,ch,volt)
      fprintf(obj.Ins,sprintf(['OUT:LIM:VOLT' num2str(ch) ' %f'],volt));
      end
      
      function LimitCurr(obj,ch,curr)
      fprintf(obj.Ins,sprintf(['OUT:LIM:VOLT' num2str(ch) ' %f'],curr));
      end
      
      function data = ReadVolt(obj,ch)
      if ch==1; fprintf(obj.Ins,'MEAS:VOLT?');
      else fprintf(obj.Ins,'MEAS:VOLT2?'); end;      
      data=str2double(fscanf(obj.Ins));
      end
            
      function data = ReadCurr(obj,ch)
      if ch==1; fprintf(obj.Ins,'MEAS:CURR?');
      else fprintf(obj.Ins,'MEAS:CURR2?'); end; 
      data=str2double(fscanf(obj.Ins));
      end
      
      function SetTimer(obj,h,m,s)
      fprintf(obj.Ins,sprintf('TIMER:HOUR %f',h));
      fprintf(obj.Ins,sprintf('TIMER:MIN %f',m));
      fprintf(obj.Ins,sprintf('TIMER:SEC %f',s));
      end
      
      function Timer(obj,sw)
      fprintf(obj.Ins,sprintf('TIMER %f',sw));
      end
   end
end