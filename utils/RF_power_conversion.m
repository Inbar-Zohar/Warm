function out = RF_power_conversion(in,direction)
% The function converts :
%   'Vpp2dBm'   Vpp <=> dBm
%   'dBm2Vpp'   dBm <=> Vpp
%   'dBm2mW'    dBm <=> mW
%   'mW2dBm'     mW <=> dBm
%   'mW2Vpp'     mW <=> Vpp
%   'Vpp2mW'    Vpp <=> mW
%   'dBV2VRMS'    dBV <=> VRMS
%   'VRMS2dBV'    VRMS <=> dBV 
if nargin<2 % default is dBm2Vpp
    direction = 'Vpp2dBm';
end

switch direction
    case 'Vpp2dBm'
        out = 10*log10(((in/2/sqrt(2)).^2)/50*1e3) ;% dBm
    case 'dBm2Vpp'
        out = 2*(2*10.^(in/10)*1e-3*50).^0.5; % Vpp
    case 'dBm2mW'
        out = 10.^(in/10); % mW
    case 'mW2dBm'
        out = 10*log10(in); % dBm
    case 'mW2Vpp'
        out = 2*(2*in*1e-3*50).^0.5 ; %Vpp
    case 'Vpp2mW'
        out = (in/2/sqrt(2)).^2/50*1e3;% mW
    case 'VRMS2dBV'
        out = 20*log10(in); % dBm
    case 'dBV2VRMS'
        out = 10.^(in/20); % dBm
end
