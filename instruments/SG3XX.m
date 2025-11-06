classdef SG3XX < handle
    %this class contains all operations for the SG38X/SG39X MW generators
    %
    
    properties
        Ins;
        Freq %Frequency in Hz
        Ampl %Amplitude in dBm
        RunStatus  %RF broadcasting?
        ModStatus  %Does the modulation work?
        ModType  %Which modulation? see Mod_type for details
        ModDeviation  %What is the modulation deviation in Hz
        ModFunction  %What is the modulation function? see Sweep for details
        ModRate  %What is the modulation rate in Hz
        ModDepth  %AM depth in percentile%
    end
    
    methods
        function obj = SG3XX(address)
            %MUST RUN AFTER CONSTRUCTOR
            %start communication and read parameters.
            obj.Ins = tcpip(address,5025,'TimeOut',2);
            fopen(obj.Ins);
            obj.getFreq; %Frequency in Hz
            obj.getAmpl; %Amplitude in dBm
            obj.RunStatus = str2double(query(obj.Ins,'ENBR?')); %RF broadcasting?
            obj.ModStatus = str2double(query(obj.Ins,'MODL?')); %Does the modulation work?
            obj.ModType = str2double(query(obj.Ins,'Type?')); %Which modulation? see Mod_type for details
            switch obj.ModType %Choose Modulation type. 0=AM ,1=FM, 2=PHASE, 3=SWEEP, 4=PULSE, 5=BLANk, 6=IQ
                case 0
                    obj.ModType = 'AM'   ; 
                case 1
                    obj.ModType = 'FM'   ; 
                case 2 
                    obj.ModType = 'PHASE';
                case 3 
                    obj.ModType = 'SWEEP'; 
                case 4 
                    obj.ModType = 'PULSE'; 
                case 5 
                    obj.ModType = 'BLANK';
                case 6 
                    obj.ModType = 'IQ'   ;
            end
            obj.ModDeviation = str2double(query(obj.Ins,'SDEV?')); %What is the modulation deviation in Hz
            obj.ModFunction = str2double(query(obj.Ins,'SFNC?')); %What is the modulation function? see Sweep for details
            switch obj.ModFunction %Function: 0=sine, 1=ramp, 2=triangle, 4=noise, 5=external
                case 0 
                    obj.ModFunction = 'SINE'    ; 
                case 1 
                    obj.ModFunction = 'RAMP'    ;
                case 2 
                    obj.ModFunction = 'TRIANGLE';
                case 4
                    obj.ModFunction = 'NOISE';
                case 5 
                    obj.ModFunction = 'EXTERNAL';
            end
            obj.ModRate = str2double(query(obj.Ins,'SRAT?')); %What is the modulation rate in Hz
            obj.ModDepth = str2double(query(obj.Ins,'ADEP?')); %What is the modulation depth in percentile%
            fclose(obj.Ins);
        end
        
        function set(obj,Freq,Ampl,unit_freq,unit_ampl)
            obj.setFreq(Freq,unit_freq)
            obj.setAmp(Ampl,unit_ampl)
        end
        
        function setFreq(obj,Freq,unit_freq)
            if ~exist('unit_freq','var'); unit_freq = 'Hz'; end
            fprintf(obj.Ins,strcat('FREQ ',num2str(Freq,'%10.10f'),' ',unit_freq));
            getFreq(obj);
        end
        
        function setAmp(obj,Ampl,unit_ampl)
            if ~exist('unit_ampl','var'); unit_freq = 'dBm'; end
            fprintf(obj.Ins,strcat('AMPR ',num2str(Ampl,'%10.10f'),' ',unit_ampl));
        end
        
        function freq = getFreq(obj)
            freq = str2double(query(obj.Ins,'FREQ?'));
            obj.Freq = freq;
        end
        
        function ampl = getAmpl(obj)
            ampl = str2double(query(obj.Ins,'AMPR?'));
            obj.Ampl = ampl;
        end

        
        function idn = IDN(obj); idn = query(obj.Ins,'*IDN?'); end
        
        function err = error(obj,display_flag)
            fprintf(obj.Ins,'LERR?');
            err = fscanf(obj.Ins);
            if exist('display_flag','var') && display_flag; disp(err); end
        end
        
        function OutputON(obj); obj.Start; end
        function OutputOFF(obj); obj.Stop; end
        
        function Start(obj)
            %Start MW broadcasting
            fprintf(obj.Ins,'ENBR 1');
            obj.RunStatus = 1;
        end
        
        function Stop(obj)
            %Stop MW broadcasting
            fprintf(obj.Ins,'ENBR 0');
            obj.RunStatus = 0;
        end
        
        function Mod_on(obj)
            %Start modulation
            fprintf(obj.Ins,'MODL 1');
            obj.ModStatus = 1;
        end
        
        function Mod_off(obj)
            %Stop modulation
            fprintf(obj.Ins,'MODL 0');
            obj.ModStatus = 0;
        end
        
        function Mod_type(obj,Mod)
            if ischar(Mod)
                switch lower(Mod); case 'am'; Mod=0; case 'fm'; Mod=1;
                    case {'pm','phase'}; Mod=2; case 'sweep'; Mod=3;
                    case 'pulse'; Mod=4; case 'blank'; Mod=5; case 'iq'; Mod=6;
                end
            end
            %Choose Modulation type. 0=AM ,1=FM, 2=PHASE, 3=SWEEP, 4=PULSE, 5=BLANk, 6=IQ
            fprintf(obj.Ins,strcat('TYPE ',int2str(Mod)));
            obj.ModType = Mod;
        end
        
        function Mod_depth(obj,d)
            %set modulation depth, d is in percetile% (0-100)
            fprintf(obj.Ins,['ADEP ' num2str(d)]);
            obj.ModDepth = d;
        end
        
        function External_IQ(obj)
            %set the IQ for external
            fprintf(obj.Ins,'QFNC 5');
            obj.ModFunction = str2double(query(obj.Ins,'SFNC?')); %What is the modulation function? see Sweep for details
            switch obj.ModFunction %Function: 0=sine, 1=ramp, 2=triangle, 5=external
                case 0; obj.ModFunction = 'SINE'    ; case 1; obj.ModFunction = 'RAMP'    ;
                case 2; obj.ModFunction = 'TRIANGLE'; case 5; obj.ModFunction = 'EXTERNAL';
            end
        end
        
        function Mod_external(obj)
            %set modulation dource as external
            fprintf(obj.Ins,'MFNC 5');
            obj.ModFunction = 'EXTERNAL';
        end
        
        function FM_Noise(obj, bw, rms)
            % Set modulation function to noise
            fprintf(obj.Ins,'MFNC 4');
            obj.ModFunction = 'NOISE';
            
            % Set noise BW (rate) and RMS (deviation)
            % NOTE: Acutal values are ?64! 
            fprintf(obj.Ins,['RATE ',num2str(bw)]);
            obj.ModRate = bw;
            fprintf(obj.Ins,['FNDV ',num2str(rms)]);
            obj.ModDeviation = rms;
        end
        
        function Sweep(obj,Deviation,Function,Rate)
            %Sweep properties
            %Function: 0=sine, 1=ramp, 2=triangle, 5=external
            %Deviation and Rate in Hz
            if ischar(Function)
                switch lower(Function)
                    case 'sine'    ; Function=0; case 'ramp'    ; Function=1;
                    case 'triangle'; Function=2; case 'external'; Function=5;
                end
            end
            fprintf(obj.Ins,strcat('SDEV ',int2str(Deviation)));
            fprintf(obj.Ins,strcat('SFNC ',int2str(Function)));
            fprintf(obj.Ins,strcat('SRAT ',int2str(Rate)));
            obj.ModDeviation = Deviation;
            obj.ModFunction = Function;
            switch obj.ModFunction %Function: 0=sine, 1=ramp, 2=triangle, 5=external
                case 0; obj.ModFunction = 'SINE'    ; case 1; obj.ModFunction = 'RAMP'    ;
                case 2; obj.ModFunction = 'TRIANGLE'; case 5; obj.ModFunction = 'EXTERNAL';
            end
            obj.ModRate = Rate;
        end
    end
    
end

