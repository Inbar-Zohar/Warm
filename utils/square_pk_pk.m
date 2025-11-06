close all;
clc;

t = 0:1e-3:12;
T = 4;
tau = T/40;

s = square(2*pi*t/T);
x = ((1 - exp(-(t-T/2*floor(2*t/T))/tau)) .* (s + 1) - exp(-(t-T/2*floor(2*t/T))/tau) .* (s - 1));
n = wgn(1, numel(t), -15);
y = x + n;
yf = medfilt1(y, 100*T);
yu = yf(yf > prctile(yf,90));
yd = yf(yf < prctile(yf,10));
Au = mode(yu);
Ad = mode(yd);

figure; hold on;
plot(t, s, 'LineWidth', 2);
plot(t, y, 'LineWidth', 2);
plot(t, yf, 'LineWidth', 2);
plot(t, Au*ones(1,numel(t)), 'k--', 'LineWidth', 2);
plot(t, Ad*ones(1,numel(t)), 'k--', 'LineWidth', 2);