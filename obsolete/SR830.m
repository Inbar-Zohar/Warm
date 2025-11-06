classdef SR830
% Updates:
% 01/04/19 - Roy:
%	1- changed "fprintf/fscanf" to "query"
%	2- edited the "path" property in "Save_Configuration"
%	3- added "Read_Configuration" function
%	4- edited "Sens_Gain" to read the sensitivity in SI units if no input is provided
% 19/11/19 - Roy: chaned Sense_Gain to SenseGain
% 18/03/21 - Roy: added auto phase

   properties
   Ins
   end
    methods
      function obj = SR830(address,is_visa)
      obj.Ins = gpib('NI', 0, address,'Timeout', 1.0);
      if exist('is_visa','var') && is_visa; obj.Ins = visa('agilent',['GPIB0::' num2str(address) '::INSTR']); end
      end
          
      function data = IDN(obj); data = query(obj.Ins,'*IDN?'); end % IDN
      function Clear(obj); fprintf(obj.Ins,'*CLS'); end % Clear
      function Reset(obj); fprintf(obj.Ins,'*RST'); end % Reset Configuration
      
      function data = Show_Channel_Screen(obj,ch) %Show channel screen (input: 1/2)
      data = query(obj.Ins,sprintf('OUTR? %d',ch));
      end
      
      function data = ShowX(obj) %Show X
      data = str2double(query(obj.Ins,'OUTP? 1'));
      end
      function data = ShowY(obj) %Show Y
      data = str2double(query(obj.Ins,'OUTP? 2'));
      end
      function data = ShowR(obj) %Show R
      data = str2double(query(obj.Ins,'OUTP? 3'));
      end
      function data = ShowTheta(obj) %Show omega
      data = str2double(query(obj.Ins,'OUTP? 4'));
      end
      function data = ShowPhase(obj) %Show phase
      data = str2double(query(obj.Ins,'PHAS?'));
      end
      function data = ShowFreq(obj) %Show frequency
      data = str2double(query(obj.Ins,'FREQ?'));
      end
      function data = ShowAmp(obj) %Show amplitude
      data = str2double(query(obj.Ins,'SLVL?'));
      end
      
      function Change_Display_Theta(obj) %change the disply in the device
      fprintf(obj.Ins,'DDEF 2,1,0');
      end
      function Change_Display_R(obj) %change the disply in the device
      fprintf(obj.Ins,'DDEF 1,1,0');
      end
      function Change_Phase_Mod_Int(obj) %Show phase mod to internal
      fprintf(obj.Ins,'FMOD 1');
      end
      function Change_Phase_Mod_Ext(obj) %Show phase mod to external (Ref In)
      fprintf(obj.Ins,'FMOD 0');
      end
      
      function SelectA(obj) %Select A
      fprintf(obj.Ins,'ISRC 0');
      end
      
      function data = SensGain(obj,num) %Change the gain
          if exist('num','var') && isnumeric(num) && (num>0 && num<1); num = 26+round(log10(num)*3); end
          if exist('num','var'); fprintf(obj.Ins,sprintf('SENS %d',num)); data = [];
          else, data = str2double(query(obj.Ins,'SENS?')); data = round(10^((mod(data,3)+1)/3))*10^(floor(data/3)-9);
          end
      end
      function data = SetPhase(obj,num) %Change the phase
          if exist('num','var') && isnumeric(num) && (num>0 && num<1); num = round(log10(num)*2)+10; end
          if exist('num','var'); fprintf(obj.Ins,sprintf('PHAS %d',num)); data = [];
          else, data = str2double(query(obj.Ins,'PHAS?')); data = round(10^(mod(data,2)/2))*10^(floor(data/2)-5);
          end
      end
      function data = SensTConst(obj,num) %Change the Time Const
          if exist('num','var') && isnumeric(num) && (num>0 && num<1); num = round(log10(num)*2)+10; end
          if exist('num','var'); fprintf(obj.Ins,sprintf('OFLT %d',num)); data = [];
          else, data = str2double(query(obj.Ins,'OFLT?')); data = round(10^(mod(data,2)/2))*10^(floor(data/2)-5);
          end
      end
      function AutoPhase(obj); fprintf(obj.Ins,'APHS'); end
      
      function config = Read_Configuration(obj) %Read configuration
      fprintf(obj.Ins,'?RST');
      config.phase = str2double(query(obj.Ins,'PHAS?'));
      config.source = str2double(query(obj.Ins,'FMOD?'));
      config.freq = str2double(query(obj.Ins,'FREQ?'));
      config.reftrig = str2double(query(obj.Ins,'RSLP?'));
      config.harmonic = str2double(query(obj.Ins,'HARM?'));
      config.sineout = str2double(query(obj.Ins,'SLVL?'));
      config.inputconf = str2double(query(obj.Ins,'ISRC?'));
      config.shildground = str2double(query(obj.Ins,'IGND?'));
      config.coupling = str2double(query(obj.Ins,'ICPL?'));
      config.notch = str2double(query(obj.Ins,'ILIN?'));
      config.sens = str2double(query(obj.Ins,'SENS?'));
      config.rmod = str2double(query(obj.Ins,'RMOD?'));
      config.timeconst = str2double(query(obj.Ins,'OFLT?'));
      config.lowpass = str2double(query(obj.Ins,'OFSL?'));
      config.sync = str2double(query(obj.Ins,'SYNC?'));
      end
      
      function Save_Configuration(obj,path, fname) %Save configuration
      fprintf(obj.Ins,'?RST');
      phase = str2double(query(obj.Ins,'PHAS?')); %#ok<NASGU>
      source = str2double(query(obj.Ins,'FMOD?')); %#ok<NASGU>
      freq = str2double(query(obj.Ins,'FREQ?')); %#ok<NASGU>
      reftrig = str2double(query(obj.Ins,'RSLP?')); %#ok<NASGU>
      harmonic = str2double(query(obj.Ins,'HARM?')); %#ok<NASGU>
      sineout = str2double(query(obj.Ins,'SLVL?')); %#ok<NASGU>
      inputconf = str2double(query(obj.Ins,'ISRC?')); %#ok<NASGU>
      shildground = str2double(query(obj.Ins,'IGND?')); %#ok<NASGU>
      coupling = str2double(query(obj.Ins,'ICPL?')); %#ok<NASGU>
      notch = str2double(query(obj.Ins,'ILIN?')); %#ok<NASGU>
      sens = str2double(query(obj.Ins,'SENS?')); %#ok<NASGU>
      rmod = str2double(query(obj.Ins,'RMOD?')); %#ok<NASGU>
      timeconst = str2double(query(obj.Ins,'OFLT?')); %#ok<NASGU>
      lowpass = str2double(query(obj.Ins,'OFSL?')); %#ok<NASGU>
      sync = str2double(query(obj.Ins,'SYNC?')); %#ok<NASGU>
      if path(end)~='\'; path(end+1) = '\'; end
      if ~exist('fname','var'); fname = 'sr_830config';  end
      if exist(path,'dir'); fname = [path  fname '.mat']; else fname = [fname '.mat']; end;
      save(fname,'phase','source','freq','reftrig','harmonic','sineout','inputconf','shildground','coupling','notch','sens','rmod','timeconst','lowpass','sync');
      clear phase source freq reftrig harmonic sineout inputconf shildground coupling notch sens rmod timeconst lowpass sync fname
      end
      
      function Load_Configuration(obj,path, fname) %Load configuration
      cd(path);
      if ~exist('fname','var'); fname = 'sr_830config';  end
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
      
      function C = readEntireConfig(obj)
          %% outputs:
        Ch2Output = str2double(query(obj.Ins,'FPOP ? 2'));
        switch Ch2Output; case 0; C.Ch2Output = 'Disp'; case 1; C.Ch2Output = 'Y'; end
        Ch2Disp = round(str2double(query(obj.Ins,'DDEF ? 2'))/10);
        switch Ch2Disp; case 0; C.Ch2Disp = 'Y'; case 1; C.Ch2Disp = 'theta'; case 2; C.Ch2Disp = 'Yn'; case 3; C.Ch2Disp = 'Aux3'; case 4; C.Ch2Disp = 'Aux4'; end
        Ch1Output = str2double(query(obj.Ins,'FPOP ? 1'));
        switch Ch1Output; case 0; C.Ch1Output = 'Disp'; case 1; C.Ch1Output = 'X'; end
        Ch1Disp = round(str2double(query(obj.Ins,'DDEF ? 1'))/10);
        switch Ch1Disp; case 0; C.Ch1Disp = 'X'; case 1; C.Ch1Disp = 'R'    ; case 2; C.Ch1Disp = 'Xn'; case 3; C.Ch1Disp = 'Aux1'; case 4; C.Ch1Disp = 'Aux3'; end
        d = sscanf(query(obj.Ins,'OEXP? 1'),'%d,%d'); C.OffsetX = d(1); C.ExpandX = 10^d(2);
        d = sscanf(query(obj.Ins,'OEXP? 2'),'%d,%d'); C.OffsetY = d(1); C.ExpandY = 10^d(2);
        d = sscanf(query(obj.Ins,'OEXP? 3'),'%d,%d'); C.OffsetR = d(1); C.ExpandR = 10^d(2);
          %% input
        ISource = str2double(query(obj.Ins,'ISRC ?'));
        switch ISource; case 0; C.ISource = 'A'; case 1; C.ISource = 'A-B'; case 2; C.ISource = 1e6; case 3; C.ISource = 100e6; end
        Grounding = str2double(query(obj.Ins,'IGND ?'));
        switch Grounding; case 0; C.Grounding = 'Float'; case 1; C.Grounding = 'Ground'; end
        Coupling = str2double(query(obj.Ins,'ICPL ?'));
        switch Coupling; case 0; C.Coupling = 'AC'; case 1; C.Coupling = 'DC'; end
        C.LineFiltDB = 10*str2double(query(obj.Ins,'ILIN ?'));
          %% refs
        RefSource = str2double(query(obj.Ins,'FMOD ?'));
        switch RefSource; case 0; C.RefSource = 'external'; case 1; C.RefSource = 'internal'; end
        C.RefFreq = str2double(query(obj.Ins,'FREQ ?'));
        C.RefPhase = str2double(query(obj.Ins,'PHAS ?'));
        RefWave = str2double(query(obj.Ins,'RSLP ?'));
        switch RefWave; case 0; C.RefWave = 'Sin'; case 1; C.RefWave = 'TTLrise'; case 2; C.RefWave = 'TTLfall'; end
        C.HarmonicDetected = str2double(query(obj.Ins,'HARM ?'));
        C.SineOutRMS = str2double(query(obj.Ins,'SLVL ?'));
          %% amplifier parameters
        Sensitivity = 10.^(8/24*(str2double(query(obj.Ins,'SENS ?'))-2) -8 );
        C.Sensitivity = round(Sensitivity*10^-floor(log10(Sensitivity)))*10^floor(log10(Sensitivity));
            C.OffsetX = C.OffsetX*C.Sensitivity;
            C.OffsetY = C.OffsetY*C.Sensitivity;
            C.OffsetR = C.OffsetR*C.Sensitivity;
        Reserve = str2double(query(obj.Ins,'RMOD ?'));
        switch Reserve; case 0; C.Reserve = 'High'; case 1; C.Reserve = 'Normal'; case 2; C.Reserve = 'Low'; end
        TConst = 10.^(9/18*str2double(query(obj.Ins,'OFLT ?')) -5 );
        C.TConst = round(TConst*10^-floor(log10(TConst)))*10^floor(log10(TConst));
        C.FilterSlope = 6*(str2double(query(obj.Ins,'OFSL ?'))+1);
        SyncFilter = str2double(query(obj.Ins,'SYNC ?'));
        switch SyncFilter; case 0; C.SyncFilter = 'OFF'; case 1; C.SyncFilter = 'ON'; end
          %% not treated (p. 10 / 1-8 in the manual):
% AUX output setup, General setup (COMM, alarms...), Data storage & transfer, interface, STATUS & errors
      end
      
    function C = writeEntireConfig(obj,C)
%           %% outputs:
%         switch C.Ch2Output; case 'Disp'; fprintf(obj.Ins,'FPOP 0 2'); case 'Y'; fprintf(obj.Ins,'FPOP 1 2'); end 
%         switch C.Ch2Disp; case 'Y'; fprintf(obj.Ins,'DDEF 0 2'); case 'theta'; fprintf(obj.Ins,'DDEF 1 2'); case 'Yn'; fprintf(obj.Ins,'DDEF 2 2'); case 'Aux3'; fprintf(obj.Ins,'DDEF 3 2'); case 'Aux4'; fprintf(obj.Ins,'DDEF 4 2'); end
%         switch C.Ch1Output; case 'Disp'; fprintf(obj.Ins,'FPOP 0 1'); case 'X'; fprintf(obj.Ins,'FPOP 1 1'); end
%         switch C.Ch1Disp; case 'X'; fprintf(obj.Ins,'DDEF 0 1'); case 'R'    ; fprintf(obj.Ins,'DDEF 1 1'); case 'Xn'; fprintf(obj.Ins,'DDEF 2 1'); case 'Aux1'; fprintf(obj.Ins,'DDEF 3 1'); case 'Aux2'; fprintf(obj.Ins,'DDEF 4 1'); end
%         fprintf(obj.Ins,['OEXP 1 ' num2str(C.OffsetX) ',' num2str(log10(C.ExpandX))]);
%         fprintf(obj.Ins,['OEXP 2 ' num2str(C.OffsetY) ',' num2str(log10(C.ExpandY))]);
%         fprintf(obj.Ins,['OEXP 3 ' num2str(C.OffsetR) ',' num2str(log10(C.ExpandR))]);
%           %% input
%         switch C.ISource; case 'A'; fprintf(obj.Ins,'ISRC 0'); case 'A-B'; fprintf(obj.Ins,'ISRC 1'); case 1e6; fprintf(obj.Ins,'ISRC 2'); case 100e6; fprintf(obj.Ins,'ISRC 3'); end
%         switch C.Grounding; case 'Float'; fprintf(obj.Ins,'IGND 0'); case 'Ground'; fprintf(obj.Ins,'IGND 1'); end
%         switch C.Coupling; case 'AC'; fprintf(obj.Ins,'ICPL 0'); case 'DC'; fprintf(obj.Ins,'ICPL 1'); end
%         fprintf(obj.Ins,['ILIN ' num2str(C.LineFiltDB/10)]);
%           %% refs
%         switch C.RefSource; case 'external'; fprintf(obj.Ins,'FMOD 0'); case 'internal'; fprintf(obj.Ins,'FMOD 1'); end
%         fprintf(obj.Ins,['FREQ ' num2str(C.RefFreq)]);
%         fprintf(obj.Ins,['PHAS ' num2str(C.RefPhase)]);
%         fprintf(obj.Ins,['RSLP ' num2str(RefWave)]);
%         switch C.RefWave; case 'Sin'; fprintf(obj.Ins,'RSLP 0'); case 'TTLrise';  fprintf(obj.Ins,'RSLP 1'); case 'TTLfall';  fprintf(obj.Ins,'RSLP 2'); end
%         fprintf(obj.Ins,['HARM ' num2str(C.HarmonicDetected)]);
%         fprintf(obj.Ins,['SLVL ' num2str(C.SineOutRMS)]);
%           %% amplifier parameters
%          = fprintf(obj.Ins,['SENS ']);
%         C.Sensitivity = round(Sensitivity*10^-floor(log10(Sensitivity)))*10^floor(log10(Sensitivity));
%             C.OffsetX = C.OffsetX*C.Sensitivity;
%             C.OffsetY = C.OffsetY*C.Sensitivity;
%             C.OffsetR = C.OffsetR*C.Sensitivity;
%         Reserve = str2double(query(obj.Ins,'RMOD ?'));
%         switch Reserve; case 0; C.Reserve = 'High'; case 1; C.Reserve = 'Normal'; case 2; C.Reserve = 'Low'; end
%         TConst = 10.^(9/18*str2double(query(obj.Ins,'OFLT ?')) -5 );
%         C.TConst = round(TConst*10^-floor(log10(TConst)))*10^floor(log10(TConst));
%         C.FilterSlope = 6*(str2double(query(obj.Ins,'OFSL ?'))+1);
%         SyncFilter = str2double(query(obj.Ins,'SYNC ?'));
%         switch SyncFilter; case 0; C.SyncFilter = 'OFF'; case 1; C.SyncFilter = 'ON'; end
%           %% not treated (p. 10 / 1-8 in the manual):
% % AUX output setup, General setup (COMM, alarms...), Data storage & transfer, interface, STATUS & errors
    end
    
   end
end