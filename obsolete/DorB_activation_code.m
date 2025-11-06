%% Dor B Activation code
instrreset; 
scope1 = DSOX3000('0x2A8D::0x1766::MY61265015'); fopen(scope1.Ins); %you need to change the (-) to (::)
scope1.IDN
%scope1 = DSOX3000('0x2A8D::0x1766::MY61265015'); fopen(scope1.Ins); %you need to change the (-) to (::)
AG1 = AG33500B('192.168.1.15','tcpip'); fopen(AG1.Ins);
AG1.IDN