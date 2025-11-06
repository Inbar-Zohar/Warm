%% reset old connections to instruments if exist
instrreset

%% establish new connections 
% Scopes
% scope1 = DSOX3000('0x2A8D::0x1766::MY61265015'); fopen(scope1.Ins); scope1.Ins.Timeout = 3;
scope2 = DSOX3000('0x2A8D::0x1766::MY56311506'); fopen(scope2.Ins); scope2.Ins.Timeout = 3;
% scope3 = DSOX3000('0x2A8D::0x1766::MY63420878'); fopen(scope3.Ins); scope3.Ins.Timeout = 3;
% warning('Scope 3, the one that was used to measure Xe phases and amplitudes are gone. Need to connect a new one.')
% scope4 = DSOX3000('0x0957::0x17A6::MY53160409'); fopen(scope4.Ins); scope4.Ins.Timeout = 3;

% Waveform Generators
AG1 = AG33500B('192.168.2.11','tcpip'); fopen(AG1.Ins); AG1.IDN;
AG2 = AG33500B('192.168.2.12','tcpip'); fopen(AG2.Ins); AG2.IDN;
AG3 = AG33500B('192.168.2.13','tcpip'); fopen(AG3.Ins); AG3.IDN;
AG4 = AG33500B('192.168.2.14','tcpip'); fopen(AG4.Ins); AG4.IDN;
AG5 = AG33500B('192.168.2.15','tcpip'); fopen(AG5.Ins); AG5.IDN;
AG6 = AG33500B('192.168.2.16','tcpip'); fopen(AG6.Ins); AG6.IDN;

% Frequency Generators for wavelength locking
SG1 = SG3XX('192.168.2.21'); fopen(SG1.Ins); SG1.IDN;
SG2 = SG3XX('192.168.2.22'); fopen(SG2.Ins); SG2.IDN;
SG3 = SG3XX('192.168.2.23'); fopen(SG3.Ins); SG3.IDN;

% Heating 
% try
%     Tabor1 = WW1071('192.168.2.24'); pause(1); fopen(Tabor1.Ins); %Tabor1.setSingleTrigOut; % 250mV p2p is for heating and creating the suitable vapor pressure of Rb
% catch E
%     display(E.message);
%     display('2nd tabor1 attempt');
%     Tabor1 = WW1071('192.168.2.24'); pause(1); fopen(Tabor1.Ins); %Tabor1.setSingleTrigOut;
% end
% Tabor1.IDN

% % Tempearature reader
% tempReader = KEITHLEY2700(6,1); fopen(tempReader.Ins); tempReader.IDN 
% fclose(tempReader.Ins)
% fclose(Tabor1.Ins)

% Current supplies
% CS1 = CS580('COM5'); fopen(CS1.Ins); CS1.IDN;
PR1 = PRO8000('COM6',1:3); fopen(PR1.Ins); PR1.IDN;

% LIAs
SR1 = SR86X('192.168.2.17'); fopen(SR1.Ins);
SR2 = SR86X('192.168.2.18'); fopen(SR2.Ins);
SR3 = SR86X('192.168.2.19'); fopen(SR3.Ins);
SR4 = SR830(1,15); fopen(SR4.Ins);

% PIDs
SIM1 = SIM900('COM1'); fopen(SIM1.Ins); SIM1.IDN;

% Laser Controllers

disp('instruments connected')


%%
DorB = struct;

DorB.Bx = GenericDevice(AG1,2);
DorB.By = GenericDevice(AG2,2); % was AG2 before changes
DorB.Bx_sim = GenericDevice(AG5,1);
DorB.By_sim = GenericDevice(AG5,2);
if exist('CS1')
    DorB.Bz = GenericDevice(CS1);
end
if exist('PR1')
    DorB.gradBx = GenericDevice(PR1,1);
    DorB.gradBy = GenericDevice(PR1,2);
    DorB.gradBz = GenericDevice(PR1,3);
end

DorB.ESRDrive = GenericDevice(AG3,1);
DorB.Xe129Drive = GenericDevice(AG1,1); % was AG1 before changes
DorB.Xe131Drive = GenericDevice(AG2,1);

if exist('SR1'); DorB.ESR_LIA = GenericDevice(SR1); end
if exist('SR2'); DorB.Xe131_LIA = GenericDevice(SR2); end
if exist('SR3'); DorB.Xe129_LIA = GenericDevice(SR3); end
if exist('SR4'); DorB.Sim_LIA = GenericDevice(SR4); end

if exist('SIM1')
    DorB.Xe129_PID = GenericDevice(SIM1,1); 
    DorB.Xe131_PID = GenericDevice(SIM1,3);     
    DorB.ESR_PID = GenericDevice(SIM1,5);  
    DorB.Bz_PID = GenericDevice(SIM1,7);  
end

DorB.Trigger = GenericDevice(AG6,2);

% DorB.PowerPump = GenericDevice(AG3,2);

DorB.LaserTempPump = GenericDevice(AG4,1);
DorB.LaserTempProbe = GenericDevice(AG4,2);

master_detuning = getMasterDetuning(87,2,1) * 1e9;
if exist('SG1')
    % DorB.DeltaPump = DetuningController(GenericDevice(SG1), 64, -1, 3.03e9);
    DorB.DeltaPump = DetuningController(GenericDevice(SG1), 64, -1, master_detuning);
end
if exist('SG2')
    % DorB.DeltaProbe = DetuningController(GenericDevice(SG2), 64, 1, 2*6e9 + 3.03e9);
    DorB.DeltaProbe = DetuningController(GenericDevice(SG2), 64, -1, master_detuning);
end
if exist('SG3')
    DorB.EOMProbe = GenericDevice(SG3);
end

%% Some defaults
DorB.ESRDrive.apply('VLim',[-0.5 0.5]);

LIAs = {'ESR_LIA','Xe129LIA','Xe131LIA'};
for ind=1:3
    LIA = LIAs{ind};
    if isfield(DorB,LIA)
        DorB.(LIA).apply('setSignalInputMode','A');
        DorB.(LIA).apply('setSignalInputCoupling','DC');
        DorB.(LIA).apply('setSignalInputShield','Ground');
        DorB.(LIA).apply('setSignalInputType','Voltage');  
        DorB.(LIA).apply('setRefSource','EXT');     
        DorB.(LIA).apply('setExtRefMode','POS');        
        DorB.(LIA).apply('setExtRefTermination','1 MOhm');
    end
end