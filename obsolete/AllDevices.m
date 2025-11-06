classdef AllDevices
    
    properties
        DeviceMap
    end
    
    methods
        function obj = AllDevices(DeviceMap)
            obj.DeviceMap = DeviceMap;
        end
        
        function varargout = apply(obj,phys_quant,func,varargin)
            param_ind = find(strcmp(obj.DeviceMap(:,1),phys_quant));
            device = obj.DeviceMap{param_ind,2}{1};
            device_type = class(device);
            nout = nargout([device_type '>' device_type '.' func]);
            varargout = cell(1,nout);

            if isa(device, 'AG33500B') || isa(device, 'PRO8000')
                % for devices with multiple channels where first input is
                % always channel/slot
                channel = obj.DeviceMap{param_ind,2}{2};
                [varargout{:}] = device.(func)(channel, varargin{:});
            elseif isa(device, 'CS580') || isa(device, 'SG3XX')
                % for devices with a single channel
                [varargout{:}] = device.(func)(varargin{:});
            else
                error(['Support for ' device_type ' not implemented yet'])
            end
        end
        
    end
end