% connect to DAQ
if ~exist('cDAQ','var')
    cDAQ = CDAQ9174;
end

% clean channels
for ind_card=1:3
    for ind_ch=0:3
        cDAQ.removeCh(ind_card, ind_ch);
    end
end
% cDAQ.removeCh(1, 0);
% cDAQ.removeCh(1, 1);
% cDAQ.removeCh(1, 2);
% cDAQ.removeCh(1, 3);
% cDAQ.removeCh(2, 0);
% cDAQ.removeCh(2, 1);
% cDAQ.removeCh(2, 2);
% cDAQ.removeCh(2, 3);
% cDAQ.removeCh(3, 0);
% cDAQ.removeCh(3, 1);
% cDAQ.removeCh(3, 2);
% cDAQ.removeCh(3, 3);

% set channels
%DEFINE THE cDAQ CHANNEL NAMES
%If required to remove channels in other scripts, just call get_daq when
%finished. Don't start to add channels and define their names in other
%places in the code. It's ugly!
cDAQ.addCh(1, 0); cDAQ.Ins.Channels(1).Name = 'PID-131 out'; 
cDAQ.addCh(1, 1); cDAQ.Ins.Channels(2).Name = 'PID-129 out'; 
cDAQ.addCh(1, 2); cDAQ.Ins.Channels(3).Name = 'PID-sim out'; 
cDAQ.addCh(1, 3); cDAQ.Ins.Channels(4).Name = 'PID-simR-out';%'Empty'; 
cDAQ.addCh(2, 0); cDAQ.Ins.Channels(5).Name = 'Xe-131 Y' ;
cDAQ.addCh(2, 1); cDAQ.Ins.Channels(6).Name = 'Xe-129 Y'; 
cDAQ.addCh(2, 2); cDAQ.Ins.Channels(7).Name = 'Xe-131 R';
cDAQ.addCh(2, 3); cDAQ.Ins.Channels(8).Name = 'Xe-129 R';
cDAQ.addCh(3, 0); cDAQ.Ins.Channels(9).Name = 'Probe Monitor PD';
cDAQ.addCh(3, 1); cDAQ.Ins.Channels(10).Name = 'Pump Monitor PD';
cDAQ.addCh(3, 2); cDAQ.Ins.Channels(11).Name = 'B sim Y';
cDAQ.addCh(3, 3); cDAQ.Ins.Channels(12).Name = 'B sim R';

cDAQ.getChannels;