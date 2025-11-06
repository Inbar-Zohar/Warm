%get noise spectrom of esr lia
%using T, ESR_LIA_out from WR mesurament. 
fs = 1 / mean(diff(T));
[psd,fff]=pwelch(ESR_LIA_out(2:end) - mean(ESR_LIA_out(2:end)), [], [], [], fs);

cd(fullfile(folder_path,'plots'));

figure ('Name', 'ESR_LIA vs time')
%using magnetometer sensitivity mag_sesns
h = plot(linspace(0,601,length(ESR_LIA_out)), ESR_LIA_out); 
xlabel('Time (sec)'); ylabel('ESR R (V)'); grid on; % time
saveas(h,'ESR_LIA_vs_time.fig')
saveas(h,'ESR_LIA_vs_time.png')

figure (); 
h = loglog(fff, sqrt(psd).*(1e-4*1e15/mag_sens)); 
grid on; xlabel('Frequency (HZ)'); ylabel('PSD of ESR LIA (fT/{\surd}Hz)'); %psd
saveas(h,'PSD_of_ESR_LIA.fig')
saveas(h,'PSD_of_ESR_LIA.png')

L=length(fff);

figure();
h = loglog(fff(L/10.3 -100:L/10.3+100), sqrt(psd(L/10.3 -100:L/10.3+100)).*(1e-4*1e15/mag_sens)); 
grid on; xlabel('Frequency (HZ)'); ylabel('PSD of ESR LIA (fT/{\surd}Hz)'); %psd around 99 Hz
saveas(h,'PSD_of_ESR_LIA_99Hz.fig')
saveas(h,'PSD_of_ESR_LIA_99Hz.png')

figure();
h = loglog(fff(L/3.0555 -100:L/3.0555+100), sqrt(psd(L/3.0555 -100:L/3.0555+100)).*(1e-4*1e15/mag_sens)); 
grid on; xlabel('Frequency (HZ)'); ylabel('PSD of ESR LIA (fT/{\surd}Hz)'); %psd around 335 Hz
saveas(h,'PSD_of_ESR_LIA_335Hz.fig')
saveas(h,'PSD_of_ESR_LIA_335Hz.png')

fft1(T,ESR_LIA_out.*(1e-4*1e15/mag_sens),0,1);%fft_of_ESR_LIA
cd('C:\Users\Constantine\SynologyDrive\Dor B\DorB_operational');
