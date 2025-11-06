T = meas.cDAQ_t;
V = meas.cDAQ_v/mean(meas.cDAQ_v);
fs = 1 / mean(diff(T));
[f,P1,pks,locks,Y] = fft1(T,V,1);
[psd,fff]=pwelch(V(2:end) - mean(V(2:end)), [], [], [], fs);
srPSD = sqrt(psd);


%time domain
myfig(6);
spl2_1 = mysubplot(3,1,1); cla;
grid on;
box on;
plot(T, V);
xlabel('Time (sec)');
ylabel('Voltage (V)');


%fft
spl2_2 = mysubplot(3,1,2); cla;
grid on;
box on;
plot(f,P1)
xlabel('Frequency (HZ)');
ylabel('FFT (V)')

%psd
spl2_3 = mysubplot(3,1,3); cla;
grid on;
box on;
loglog(fff,srPSD)
xlabel('Frequency (HZ)');
ylabel('srPSD of (V/{\surd}Hz)')