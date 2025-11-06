classdef SR844
   properties
   Ins
   end
    methods
      function obj = SR844(address,GPIB_address)
      if exist('GPIB_address','var') && GPIB_address==0
      obj.Ins = gpib('NI', 0, address,'Timeout', 1.0);
      else
      obj.Ins = gpib('agilent', 7, address,'Timeout', 1.0);
      end
      end
      
          
      function data = IDN(obj) %IDN
          fprintf(obj.Ins,'*IDN?');
          data = fscanf(obj.Ins);
      end
      function data = Show_Channel_Screen(obj,ch) %Show channel screen (input: 1/2)
      fprintf(obj.Ins,sprintf('OUTR? %d',ch));
      data = fscanf(obj.Ins);
      end
      function Change_Display_R(obj) %change the disply in the device
      fprintf(obj.Ins,'DDEF 1,1,0');
      end
      function Change_Display_Theta(obj) %change the disply in the device
      fprintf(obj.Ins,'DDEF 2,1,0');
      end
      function Change_Output_R(obj) %change the *output* in the device
      fprintf(obj.Ins,'FPOP 1,0');
      fprintf(obj.Ins,'DDEF 1,1');
      end
      function Change_Output_X(obj) %change the *output* in the device
      fprintf(obj.Ins,'FPOP 1,1');
      end
      function Change_Output_Theta(obj) %change the *output* in the device
      fprintf(obj.Ins,'FPOP 2,0');
      fprintf(obj.Ins,'DDEF 2,1');
      end
      function Change_Output_Y(obj) %change the *output* in the device
      fprintf(obj.Ins,'FPOP 2,1');
      end
      function data = ShowX(obj) %Show X
      fprintf(obj.Ins,'OUTP? 1');
      data = str2double(fscanf(obj.Ins));
      end
      function data = ShowY(obj) %Show Y
      fprintf(obj.Ins,'OUTP? 2');
      data = str2double(fscanf(obj.Ins));
      end
      function data = ShowR(obj) %Show R
      fprintf(obj.Ins,'OUTP? 3');
      data = str2double(fscanf(obj.Ins));
      end
      function data = ShowTheta(obj) %Show omega
      fprintf(obj.Ins,'OUTP? 4');
      data = str2double(fscanf(obj.Ins));
      end
      function data = ShowPhase(obj) %Show phase
      fprintf(obj.Ins,'PHAS?');
      data = str2double(fscanf(obj.Ins));
      end
      function data = ShowFreq(obj) %Show frequency
      fprintf(obj.Ins,'FREQ?');
      data = str2double(fscanf(obj.Ins));
      end
      function data = ShowAmp(obj) %Show amplitude
      fprintf(obj.Ins,'SLVL?');
      data = str2double(fscanf(obj.Ins));
      end
      function Change_Phase_Mod_Int(obj) %Show phase mod to internal
      fprintf(obj.Ins,'FMOD 1');
      end
      function Change_Phase_Mod_Ext(obj) %Show phase mod to external (Ref In)
      fprintf(obj.Ins,'FMOD 0');
      end
      
      function Change_Phase(obj,phase)
      fprintf(obj.Ins,['PHAS',num2str(phase)]);
      end
      
      function time = elapsedTime(obj)
      if nargout<1; fprintf(obj.Ins,'SETL');
      else time = obj.SensTConst*num2str(query(obj.Ins,'SETL ?'));
      end
      end
      
      function SelectA(obj) %Select A
      fprintf(obj.Ins,'ISRC 0');
      end
      function Clear(obj) %Clear
      fprintf(obj.Ins,'*CLS');
      end
      function Reset(obj) %Reset Configuration
      fprintf(obj.Ins,'*RST');
      end
      function Sens_Gain(obj,num) %Change the gain
      fprintf(obj.Ins,sprintf('SENS %d',num));
      end
      function Save_Configuration(obj,path,fname) %Save configuration
      fprintf(obj.Ins,'?RST');
      fprintf(obj.Ins,'PHAS?');
      phase = str2double(fscanf(obj.Ins));
      fprintf(obj.Ins,'FMOD?');
      source = str2double(fscanf(obj.Ins));
      fprintf(obj.Ins,'FREQ?');
      freq = str2double(fscanf(obj.Ins));
      fprintf(obj.Ins,'RSLP?');
      reftrig = str2double(fscanf(obj.Ins));
      fprintf(obj.Ins,'HARM?');
      harmonic = str2double(fscanf(obj.Ins));
      fprintf(obj.Ins,'SLVL?');
      sineout = str2double(fscanf(obj.Ins));
      fprintf(obj.Ins,'ISRC?');
      inputconf = str2double(fscanf(obj.Ins));
      fprintf(obj.Ins,'IGND?');
      shildground = str2double(fscanf(obj.Ins));
      fprintf(obj.Ins,'ICPL?');
      coupling = str2double(fscanf(obj.Ins));
      fprintf(obj.Ins,'ILIN?');
      notch = str2double(fscanf(obj.Ins));
      fprintf(obj.Ins,'SENS?');
      sens = str2double(fscanf(obj.Ins));
      fprintf(obj.Ins,'RMOD?');
      rmod = str2double(fscanf(obj.Ins));
      fprintf(obj.Ins,'OFLT?');
      timeconst = str2double(fscanf(obj.Ins));
      fprintf(obj.Ins,'OFSL?');
      lowpass = str2double(fscanf(obj.Ins));
      fprintf(obj.Ins,'SYNC?');
      sync = str2double(fscanf(obj.Ins));
      cd(path);
      if ~exist('fname','var'); fname = 'sr_844config';  end
      save([fname '.mat'],'phase','source','freq','reftrig','harmonic','sineout','inputconf','shildground','coupling','notch','sens','rmod','timeconst','lowpass','sync');
      clear phase source freq reftrig harmonic sineout inputconf shildground coupling notch sens rmod timeconst lowpass sync fname
      end
      
      function Load_Configuration(obj,path, fname) %Load configuration
          cd(path);
          if ~exist('fname','var'); fname = 'sr_844config';  end
          S = load([fname '.mat']);
          fprintf(obj.Ins,'?RST');
          fprintf(obj.Ins,sprintf('PHAS %f',S.phase));
          fprintf(obj.Ins,sprintf('FMOD %f',S.source));
          fprintf(obj.Ins,sprintf('FREQ %f',S.freq));
          fprintf(obj.Ins,sprintf('RSLP %f',S.reftrig));
          fprintf(obj.Ins,sprintf('HARM %f',S.harmonic));
          fprintf(obj.Ins,sprintf('SLVL %f',S.sineout));
          fprintf(obj.Ins,sprintf('ISRC %f',S.inputconf));
          fprintf(obj.Ins,sprintf('IGND %f',S.shildground));
          fprintf(obj.Ins,sprintf('ICPL %f',S.coupling));
          fprintf(obj.Ins,sprintf('ILIN %f',S.notch));
          fprintf(obj.Ins,sprintf('SENS %f',S.sens));
          fprintf(obj.Ins,sprintf('RMOD %f',S.rmod));
          fprintf(obj.Ins,sprintf('OFLT %f',S.timeconst));
          fprintf(obj.Ins,sprintf('OFSL %f',S.lowpass));
          fprintf(obj.Ins,sprintf('SYNC %f',S.sync));
          clear S fname
      end
      
      function data = SensGain(obj,num) %Change the gain
          if exist('num','var') && isnumeric(num) && (num>0 && num<1); num = 26+round(log10(num)*3); end;
          if exist('num','var'); fprintf(obj.Ins,sprintf('SENS %d',num)); data = [];
          else data = str2double(query(obj.Ins,'SENS?')); data = round(10^((mod(data,3)+1)/3))*10^(floor(data/3)-9);
          end
      end
      function data = SensTConst(obj,num) %Change the gain
          if exist('num','var') && isnumeric(num) && (num>0 && num<1); num = round(log10(num)*2)+10; end;
          if exist('num','var'); fprintf(obj.Ins,sprintf('OFLT %d',num)); data = [];
          else
              data = str2double(query(obj.Ins,'OFLT?'));
          data = round(10^(mod(data,2)/2))*10^(floor(data/2)-5);
          end
      end
      
      function data = SensTslope(obj) %Change the gain
%           if exist('num','var') && isnumeric(num) && (num>0 && num<1); num = round(log10(num)*2)+10; end;
%           if exist('num','var'); fprintf(obj.Ins,sprintf('OFLT %d',num)); data = [];
%           else
              data = str2double(query(obj.Ins,'OFSL?'));
%           data = round(10^(mod(data,2)/2))*10^(floor(data/2)-5);
%           end
      end
      function C = readEntireConfig(obj)
          %% outputs:
        Ch2Output = str2double(query(obj.Ins,'FPOP ? 2'));
        switch Ch2Output; case 0; C.Ch2Output = 'Disp'; case 1; C.Ch2Output = 'Y'; end 
        Ch2Disp = round(str2double(query(obj.Ins,'DDEF ? 2'))/10);
        switch Ch2Disp; case 0; C.Ch2Disp = 'Y'; case 1; C.Ch2Disp = 'theta'; case 2; C.Ch2Disp = 'Yn'; case 3; C.Ch2Disp = 'YndBm'; case 4; C.Ch2Disp = 'Aux2'; end
        Ch1Output = str2double(query(obj.Ins,'FPOP ? 1'));
        switch Ch1Output; case 0; C.Ch1Output = 'Disp'; case 1; C.Ch1Output = 'X'; end 
        Ch1Disp = round(str2double(query(obj.Ins,'DDEF ? 1'))/10);
        switch Ch1Disp; case 0; C.Ch1Disp = 'X'; case 1; C.Ch1Disp = 'R'    ; case 2; C.Ch1Disp = 'RdBm'; case 3; C.Ch1Disp = 'Xn'; case 4; C.Ch1Disp = 'Aux1'; end
        ratio = str2double(query(obj.Ins,'DRAT ?'));
        switch ratio; case 0; C.ratio = 'none'; case 1; C.ratio = 'Aux1'; case 2; C.ratio = 'Aux2'; end
        
        C.OffsetX = str2double(query(obj.Ins,'DOFF? 1,0'));
        C.OffsetR = str2double(query(obj.Ins,'DOFF? 1,1'));
        C.OffsetRdBm = str2double(query(obj.Ins,'DOFF? 1,2'))*200;
        C.OffsetY = str2double(query(obj.Ins,'DOFF? 2,0'));
        C.ExpandX = 10^str2double(query(obj.Ins,'DEXP? 1,0'));
        C.ExpandR = 10^str2double(query(obj.Ins,'DOFF? 1,1'));
        C.ExpandRdBm = 10^str2double(query(obj.Ins,'DOFF? 1,2'));
        C.ExpandXn = 10^str2double(query(obj.Ins,'DEXP? 1,3'));
        C.ExpandY = 10^str2double(query(obj.Ins,'DEXP? 2,0'));
        C.ExpandTheta = 10^str2double(query(obj.Ins,'DEXP? 2,1'));
        C.ExpandYn = 10^str2double(query(obj.Ins,'DEXP? 2,2'));
        
          %% input
        InImpedance = str2double(query(obj.Ins,'INPZ ?'));
        switch InImpedance ; case 0; C.InImpedance = 50 ; case 1; C.InImpedance = 10e3; end
        
        C.AuxIn1Volts  = str2double(query(obj.Ins,'AUXI ? 1'));
        C.AuxIn2Volts  = str2double(query(obj.Ins,'AUXI ? 2'));
        C.AuxOut1Volts = str2double(query(obj.Ins,'AUXO ? 1'));
        C.AuxOut2Volts = str2double(query(obj.Ins,'AUXO ? 2'));
          %% refs
        RefSource = str2double(query(obj.Ins,'FMOD ?'));
        switch RefSource; case 0; C.RefSource = 'external'; case 1; C.RefSource = 'internal'; end
        C.RefFreq = str2double(query(obj.Ins,'FREQ ?'));
        C.RefFrAq = str2double(query(obj.Ins,'FRAQ ?'));
        C.IF_Freq = str2double(query(obj.Ins,'FRIQ ?'));
        C.RefPhase = str2double(query(obj.Ins,'PHAS ?'));
        C.HarmonicDetected = str2double(query(obj.Ins,'HARM ?'));
        RefImpedance = str2double(query(obj.Ins,'REFZ ?'));
        switch RefImpedance; case 0; C.RefImpedance = 50; case 1; C.RefImpedance = 1e6; end
          %% amplifier parameters
        Sensitivity = 10.^(3/6*(str2double(query(obj.Ins,'SENS ?'))-14) );
        C.Sensitivity = round(Sensitivity*10^-floor(log10(Sensitivity)))*10^floor(log10(Sensitivity));
            C.OffsetX = C.OffsetX*C.Sensitivity;
            C.OffsetR = C.OffsetR*C.Sensitivity;
            C.OffsetY = C.OffsetY*C.Sensitivity;
        WideReserve   = str2double(query(obj.Ins,'WRSV ?'));
        switch WideReserve  ; case 0; C.WideReserve   = 'High'; case 1; C.WideReserve   = 'Normal'; case 2; C.WideReserve   = 'Low'; end
        ClosedReserve = str2double(query(obj.Ins,'CRSV ?'));
        switch ClosedReserve; case 0; C.ClosedReserve = 'High'; case 1; C.ClosedReserve = 'Normal'; case 2; C.ClosedReserve = 'Low'; end
        TConst = 10.^(3/6*(str2double(query(obj.Ins,'OFLT ?'))-8) );
        C.TConst = round(TConst*10^-floor(log10(TConst)))*10^floor(log10(TConst));
        C.FilterSlope = 6*(str2double(query(obj.Ins,'OFSL ?'))+1);
        % SETL ? = measure internal timer (see elapsedTime)
          %% not treated (p. 10 / 1-8 in the manual):
% AUX output setup, General setup (COMM, alarms...), Data storage & transfer, interface, STATUS & errors
      end
    
    end
end