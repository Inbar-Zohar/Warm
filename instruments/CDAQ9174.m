classdef CDAQ9174 < dynamicprops
    % Updates:
    % 08/03/20 - ROY: Created
    
    properties(Transient)
        Ins;
        ListenerHandle;
    end
    properties
        storage
    end
    methods
        function obj = CDAQ9174
            obj.Ins = daq.createSession('ni');
            for card = 1:1; obj.addCh(card,0:3); end
            obj.setFs(2048); % in samples/sec
            obj.setDuration(10); % in sec
            % obj.Ins.IsContinuous = true;
        end
        
        function data = getChannels(obj); data = obj.Ins.Channels; end %IDN
        function Reset(obj); fprintf(obj.Ins,'*RST'); end % Reset
        
        function addCh(obj,card,ch)
            if length(card)~=length(ch)
                for cd = card
                    for Ind = 1:length(ch)
                        wrnng = warning;
                        warning off % to supress Warning: Added channel does not support on-demand operations: only clocked operations are allowed. 
                        addAnalogInputChannel(obj.Ins,['cDAQ1Mod' num2str(cd)], ch(Ind), 'Voltage');
                        warning(wrnng)
                    end
                end
            else 
                for Ind = 1:length(ch)
                    wrnng = warning;
                    warning off % to supress Warning: Added channel does not support on-demand operations: only clocked operations are allowed.
                    addAnalogInputChannel(obj.Ins,['cDAQ1Mod' num2str(card(Ind))], ch(Ind), 'Voltage');
                    warning(wrnng)
                end
            end
        end
        
        function removeCh(obj,card,chnum)
            if exist('card','var')
                C = obj.Ins.Channels;
                if length(card)==1; card = card*ones(size(chnum)); end
                for Ind = 1:length(chnum)
                    devicestr = ['cDAQ1Mod' num2str(card(Ind))]; chstr = ['ai' num2str(chnum(Ind))];
                    for Jnd = 1:length(C)
                        if strcmp(C(Jnd).Device.ID,devicestr) && strcmp(C(Jnd).ID,chstr)
                            obj.Ins.removeChannel(Jnd);
                            C = obj.Ins.Channels;
                            break
                        end
                    end
                end
            else
                for Ind = length(chnum):-1:1
                    obj.Ins.removeChannel(Ind);
                end
            end
        end
        
        function setCh(obj,card,ch)
            C = obj.Ins.Channels;
            obj.removeCh(1:length(C));
            obj.addCh(card,ch);
        end
        
        function setFs(obj,Fs)
            obj.Ins.Rate = Fs;
        end
        %         function v = readNow(obj); v = obj.Ins.inputSingleScan; end;
        
        function setDuration(obj,T); obj.Ins.DurationInSeconds = T; end
        function ListenerHandle = measureToFile(obj,fname,file_duration_min,flag_add_timestamp)
            if nargin < 4
                flag_add_timestamp = true;
            end
            %             obj.Ins.IsContinuous = 1;
            if nargout~=1
                error(['ERROR! this function creates a listener. '...
                    'The listener MUST be saved for later deletion after '...
                    'the run is finished! Keep the handle and delete when finished']);
            end
            if ~exist('file_duration_min','var'); file_duration_min = 1; end
            if  obj.Ins.NumberOfScans<obj.Ins.Rate*60*file_duration_min
                obj.Ins.NumberOfScans=obj.Ins.Rate*60*file_duration_min;
            end
            obj.Ins.NotifyWhenDataAvailableExceeds = round(obj.Ins.Rate*60*file_duration_min);
            ListenerHandle = obj.Ins.addlistener('DataAvailable',@(src,event) ...
                obj.saveToFile(event.TimeStamps,event.Data,fname,flag_add_timestamp));
            obj.Ins.startBackground;
        end
        function ListenerHandle = prepareMeasureToFile(obj,fname,file_duration_sec,total_duration_sec)
            if nargout<1
                error(['ERROR! this function creates a listener. '...
                    'The listener MUST be saved for later deletion after '...
                    'the run is finished! Keep the handle and delete when finished']);
            end
            if ~exist('file_duration_min','var'); file_duration_sec = 1; end
            if ~exist('total_duration_min','var'); total_duration_sec = file_duration_sec; end
            if isnumeric(total_duration_sec) && ~isnan(total_duration_sec) && ~isinf(total_duration_sec)
                obj.Ins.NumberOfScans=ceil(obj.Ins.Rate*total_duration_sec);
            else, obj.Ins.IsContinuous = true;
            end
            obj.Ins.NotifyWhenDataAvailableExceeds = ceil(obj.Ins.Rate*file_duration_sec);
            %             obj.ListenerHandle = obj.Ins.addlistener('DataAvailable',@(src,event) obj.saveToFile(event.TimeStamps,event.Data,fname));
            ListenerHandle = obj.Ins.addlistener('DataAvailable',@(src,event) obj.saveToFile(event.TimeStamps,event.Data,fname));
        end
        function applyMeasureToFile(obj)
            obj.Ins.startBackground;
        end
        function ListenerHandle = backgroundLongRecording(obj,fname,file_duration_min,total_duration_min)
            if varargout~=1
                error(['ERROR! this function creates a listener. '...
                    'The listener MUST be saved for later deletion after '...
                    'the run is finished! Keep the handle and delete when finished']);
            end
            if ~exist('file_duration_min','var'); file_duration_min = 1; end
            if exist('total_duration_min','var') && isnumeric(total_duration_min) && ~isnan(total_duration_min) && ~isinf(total_duration_min)
                obj.Ins.NumberOfScans=ceil(obj.Ins.Rate*60*total_duration_min);
            elseif  obj.Ins.NumberOfScans<obj.Ins.Rate*60*file_duration_min
                obj.Ins.NumberOfScans=obj.Ins.Rate*60*file_duration_min;
            end
            obj.Ins.NotifyWhenDataAvailableExceeds = ceil(obj.Ins.Rate*60*file_duration_min);
            ListenerHandle = obj.Ins.addlistener('DataAvailable',@(src,event) ...
                obj.saveToFile(event.TimeStamps,event.Data,fname));
            obj.Ins.startBackground;
        end
        
        
        function [t,v] = measureForeground(obj); [v,t] = obj.Ins.startForeground; end
        function Stop(obj); obj.Ins.stop; end
        function Pause(obj,T); if exist('T','var'); obj.Ins.wait(T); else, obj.Ins.wait; end; end
        %         function stt = dataAvailable
        
        function saveToFile(obj,t,v,fname,flag_add_timestamp)
            if nargin < 5
                flag_add_timestamp = true;
            end
            %             [data, timestamps, ~] = read(src, src.ScansAvailableFcnCount, 'OutputFormat', 'Matrix');
            if flag_add_timestamp
                clkstr = clock;
                clkstr = [...
                    '_' num2str((clkstr(1)-2000)*1e4+clkstr(2)*1e2+clkstr(3)) ...
                    '_' num2str(clkstr(4)*1e4+clkstr(5)*1e2+round(clkstr(6)),'%06.0f') ];
                if strcmp(fname(end-3:end),'.mat'); fname = [fname(1:end-4) clkstr '.mat'];
                else, fname = [fname clkstr '.mat'];
                end
            end
            save(fname,'t','v');
            disp(['saved: ' fname]);
        end
        
        function setACcoupling(obj,channels)
            if ( ischar(channels) && strcmpi(channels,'all') ) || ...
                    ( isnumeric(channels) && (isinf(channels) || channels<0) )
                channels = 1:length(obj.Ins.Channels);
            end
            for ch = channels; obj.Ins.Channels(ch).Coupling = 'AC'; end
        end
        function setDCcoupling(obj,channels)
            if ( ischar(channels) && strcmpi(channels,'all') ) || ...
                    ( isnumeric(channels) && (isinf(channels) || channels<0) )
                channels = 1:length(obj.Ins.Channels);
            end
            for ch = channels; obj.Ins.Channels(ch).Coupling = 'DC'; end
        end
        function C = getCoupling(obj,channels)
            if ( ischar(channels) && strcmpi(channels,'all') ) || ...
                    ( isnumeric(channels) && (isinf(channels) || channels<0) )
                channels = 1:length(cdaq.Ins.Channels);
            end
            C = cell(length(channels),1); for ch = channels; C{channels==ch} = obj.Ins.Channels(ch).Coupling; end
        end
        
        function clearListener(obj,lh); delete(lh); end %#ok<INUSL>
    end
end
