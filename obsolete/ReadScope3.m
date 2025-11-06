%% data bank
% D:\NMRG\NMRG_code_and_data_history 
%% BW results
% D:\NMRG\NMRGdata\BW
%%
Freqs=[ 35 32 29 26 23 20  18 16 14 12 10 9 8 7  6  5  4 3 2 1 0.5 0.25  0.1 0.05 0.025 ];


%waitforbuttonpress  


[~,ModFreqDatat] = AG3.Output(2,[],true); 
ModFreqDatatFin = ModFreqDatat{1};

FreqItV=ModFreqDatatFin.p(1);
AmpItV=ModFreqDatatFin.p(2);

% read signal from scope
[t_, v_] = scope3.Read(channels);

figure(floor(FreqItV*1000))
hold on
plot(t_,v_)
xlabel('Time [sec]')
ylabel('Scope signal [V]')
title(['Modulation Frequency ', num2str(FreqItV), '[Hz]'])

cd D:\NMRG\NMRGdata\BW\RunBWdata1
DatName=['DatMod',num2str(FreqItV*1000),'mHz']
save(DatName,'t_','v_','ModFreqDatatFin','FreqItV')
cd C:\Users\ASAFMAN\Dropbox\genb\minime_operational

