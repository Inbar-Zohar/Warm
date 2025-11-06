% ========================================================================
%   This script isn't measuring anything, it just sets the system to the
%   ESR WP optimization
% ========================================================================

if isfield(DorB,'Xe129_PID'); DorB.Xe129_PID.apply('setOutputModeMAN'); end
if isfield(DorB,'Xe131_PID'); DorB.Xe131_PID.apply('setOutputModeMAN'); end

DorB.Xe129Drive.apply('OutputOFF');
DorB.Xe131Drive.apply('OutputOFF');

DorB.Bx.apply('Sin',177,0.1,0,prm.WP.fields.Bx.B0);
DorB.Bx.apply('Syncio',1);

DorB.By.apply('Sin',220,0.1,0,prm.WP.fields.By.B0);
DorB.By.apply('Syncio',1);

DorB.ESRDrive.apply('OutputON');