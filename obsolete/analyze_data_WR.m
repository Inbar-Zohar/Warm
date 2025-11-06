%% Load data
clear
load('D:\Dor_B\Dor_B_code_and_data_history\2023\12\18\2023_12_18__12_46_55\experiment_raw_data.mat')
%% Plot 

smooth_factor = 1;
t = meas.tp_wrld_rot;
Fs = 1/(t(2)-t(1));

temp = meas.vp_wrld_rot;
i = 11; 
probe_voltage = smooth(temp(i,:),smooth_factor);
i = 3; 
Xe129_R = smooth(temp(i,:),smooth_factor);
i = 4; 
Xe131_R = smooth(temp(i,:),smooth_factor);

figure(1), clf
hold on
plot(t,normalize(probe_voltage))
plot(t,normalize(Xe129_R))
plot(t,normalize(Xe131_R))
hold off
xlabel('Time [sec]')


[xcorr_a_b,lag]= xcorr(normalize(Xe129_R),normalize(Xe131_R));
[~,I] = max(abs(xcorr_a_b));
lagDiff = lag(I);
timeDiff = lagDiff/Fs;

[xcorr_a_a,lag]= xcorr(normalize(Xe129_R),normalize(Xe129_R));
[xcorr_b_b,lag]= xcorr(normalize(Xe131_R),normalize(Xe131_R));


figure(2); clf; hold all; %
plot(lag/Fs,xcorr_a_b/max(xcorr_a_b))
plot(lag/Fs,xcorr_a_a/max(xcorr_a_a))
plot(lag/Fs,xcorr_b_b/max(xcorr_b_b))

function out = normalize(x)

out = x-mean(x(:));
out = out/max(abs(out(:)));
% out=x;
end