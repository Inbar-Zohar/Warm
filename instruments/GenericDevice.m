classdef GenericDevice
    % This is a wrapper class for a generic single output. It recieved the 
    % device object and, if relevant, the channel. This can be used to assign
    % a single output to the physical quantity it controls. 
    % It can call in function of the object (without need to specify
    % channel) by calling the "apply" command, with the first argument the
    % name of the desired function as a string and followed by the
    % arguments of the function as usual. 

    properties
        device
        channel
    end
    
    methods
        function obj = GenericDevice(device,channel)
            obj.device = device;

            if isa(device, 'AG33500B') || isa(device, 'PRO8000') || isa(device,'SIM900')
                if nargin < 2
                    error(['Need to define a channel for device_type ' obj.device_type])
                else
                    obj.channel = channel;
                end
            else
                obj.channel = nan;
            end
        end
        
        function varargout = apply(obj,func,varargin)
            device_type = class(obj.device);

            if ~ismethod(obj.device,func); error("Function %s doesn't exist for device %s", func, device_type); end

            nout = nargout([device_type '>' device_type '.' func]);
            varargout = cell(1,nout);

            if isfinite(obj.channel)
                % for devices with multiple channels where first input is
                % always channel/slot
                [varargout{:}] = obj.device.(func)(obj.channel, varargin{:});
            else
                % for devices with a single channel
                [varargout{:}] = obj.device.(func)(varargin{:});
            end
        end
        
    end
end