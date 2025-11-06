function [f, P1, pks, locs, Y] = fft1(t,v,get_peaks,to_plot)

% check for input parameters
switch nargin
    case 2
        get_peaks = 0;      % return peaks
        to_plot   = 0;      % plot the single-sided FFT
    
    case 3
        to_plot   = 0;
%         get_peaks = 1;
    
    case 4
%         to_plot   = 1;
%         get_peaks = 1;
end


% t = (0:L-1)*T;        % Time vector
L = length(v);%1000;             % Length of signal
T = mean(diff(t));%1/Fs;             % Sampling period
Fs = 1/T;            % Sampling frequency
S = v-mean(v);%0.7*sin(2*pi*50*t) + sin(2*pi*120*t);

X = S;% + 2*randn(size(t));
% X = v;
% figure(1)
% plot(t,X)
% title('Signal Corrupted with Zero-Mean Random Noise')
% xlabel('t (sec)')
% ylabel('X(t)')

Y = fft(X);
P2 = abs(Y/L);
P1 = P2(1:floor(L/2)+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:floor(L/2))/L;

if to_plot
    figure
    hold all
    plot(f,P1)
    title('Single-Sided Amplitude Spectrum of X(t)')
    xlabel('f (Hz)')
    ylabel('|P1(f)|')
end

pks = 0; locs = 0;
if get_peaks
%     [pks, locs, ~, ~] = findpeaks(P1,f,'MinPeakHeight', 0.1);
    [pks, locs, ~, ~] = findpeaks(P1,f,'MinPeakProminence', 0.02*max(P1));
    if to_plot
        plot(locs,pks,'o');
end

end