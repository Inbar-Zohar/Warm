classdef DetuningController 
    % This is a wrapper class for control over laser lock via a Vescent
    % OPLL D-135 box. It receives as arguments the frequency source (as a
    % GenericDevice object), the manually-set OPLL parameters םf
    % multiplier and sign, and optionally a frequency offset to correct
    % between the actual beat frequency at lock to detuning relative to the
    % relevant atomic line. 
    % 
    % Example: let's say we want to lock a laser 3 GHz to the red from the 
    % F=3 line of Rb85, and our master if locked to the F=2 line. The OPLL 
    % sign should be (-) (so the argument OPLL_sign is -1), and the
    % frequency offset is simply the hyperfine detuning (so freqOffset is 
    % 3.03e9). Therefore, when we set detuning of -3 GHz, the frequency of
    % the beatnote should be -3-3.03=-6.03 GHz, and the frequency generator 
    % should be set to -6.03 GHz / (-1 * OPLL_N).
    
    properties
        detuningDevice
        OPLL_N
        OPLL_sign
        freqOffset
        minFreq
        maxFreq
    end
    
    methods
        function obj = DetuningController(detuningDevice,OPLL_N,OPLL_sign,freqOffset)
            obj.detuningDevice = detuningDevice;
            obj.OPLL_N = OPLL_N;
            obj.OPLL_sign = OPLL_sign;
            if nargin < 4; freqOffset = 0; end
            obj.freqOffset = freqOffset;
            obj.minFreq = 30e6;
            obj.maxFreq = 240e6;            
        end

        function setDetuning(obj, detuning)
            freq_to_set = (detuning - obj.freqOffset) / (obj.OPLL_sign * obj.OPLL_N);
            if freq_to_set >= obj.minFreq && freq_to_set <= obj.maxFreq
                if isa(obj.detuningDevice.device,'SG3XX')
                    obj.detuningDevice.apply('setFreq',freq_to_set);
                else
                    error(['Support for detuning device ' class(obj.detuningDevice.device) ...
                        ' not implemeneted yet']);
                end
            else
                error(['Requested detuning of ' num2str(detuning/1e6) ...
                    ' MHz requires frequency of ' num2str(freq_to_set/1e6) ...
                    ' MHz which is outside of the allowed range of ' ...
                    num2str(obj.minFreq/1e6) '-' num2str(obj.maxFreq/1e6) ' MHz.'])
            end
        end
        
        function detuning = getDetuning(obj)
            if isa(obj.detuningDevice.device,'SG3XX')
                freq = obj.detuningDevice.apply('getFreq');
            else
                error(['Support for detuning device ' class(obj.detuningDevice.device) ...
                    ' not implemeneted yet']);
            end
            detuning = freq * obj.OPLL_N * obj.OPLL_sign + obj.freqOffset;
        end        
    end
end