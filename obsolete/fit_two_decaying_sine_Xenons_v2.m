function [cf_,gof] = fit_two_decaying_sine_Xenons_v2(t, v, is_fit_cos, majority_flag, to_plot)
% Edited at 31/10/2021 by Roy Elkabetz 
% set variables
if ~exist('majority_flag', 'var') || ~isnumeric(majority_flag) || any(majority_flag==[129, 131])
    majority_flag = 129;
end
if size(t, 1)~=1 
    t = t.'; 
end 
if size(v, 1)~=1 
    v = v.'; 
end
if ~exist('to_plot', 'var') 
    to_plot = true; 
end

% set signal indices limits
i_min = 1;
i_max = length(t);

% set length of fft padding
N_padd = 2^21;

BW_POINTS = 5;

% subtract Dc from signal
sig_fft = v(i_min:i_max) - mean(v(i_min:i_max));

% compute the n point DFT of the processed signal
yyy = abs(fftshift(fft(sig_fft, N_padd)));

% get the angular frequency vector
w_vec = 2 * pi * linspace(-0.5 / (t(2) - t(1)), 0.5 / (t(2) - t(1)), N_padd);

mask_yyy = ones(1, N_padd); 
mask_yyy(round(N_padd / 2) - BW_POINTS:round(N_padd / 2) + BW_POINTS) = 0;
mask_yyy(1:100) = 0; 
mask_yyy(end - 100:end) = 0; 
mask_yyy(w_vec <= 0) = 0;
ind_omega = find(mask_yyy .* yyy == max(mask_yyy .* yyy), 1);
omega1 = abs(w_vec(ind_omega));

x_approx = 3.37;

omega2A = omega1 / x_approx;
mask2A = abs(w_vec / omega2A - 1) < 0.5;
[A2A, ind2A] = max(mask2A .* yyy);
omega2Af = abs(w_vec(ind2A));

omega2B = omega1 * x_approx;
mask2B = abs(w_vec / omega2B - 1) < 0.5;
[A2B, ind2B] = max(mask2B .* yyy);
omega2Bf = abs(w_vec(ind2B));

if A2A > A2B && abs(omega2Af / omega2A - 1) < 0.9 / 2
    omega131g = omega2Af; 
    Ind131 = ind2A; 
    omega129g = omega1; 
    Ind129 = ind_omega;
elseif A2B > A2A && abs(omega2Bf / omega2B - 1) < 0.9 / 2
    omega129g = omega2Bf; 
    Ind129 = ind2B; 
    omega131g = omega1; 
    Ind131 = ind_omega;
elseif majority_flag == 129
    omega129g = omega1; 
    Ind129 = ind_omega;
    omega131g = omega129g / x_approx;
    Ind131 = find(abs(w_vec - omega131g) == min(abs(w_vec - omega131g)));
elseif majority_flag == 131
    omega131g = omega1; 
    Ind131 = ind_omega;
    omega129g = omega131g * x_approx;
    Ind129 = find(abs(w_vec - omega129g) == min(abs(w_vec - omega129)));
else; error('no two distinct frequency guesses');
end

FWHM_max = 1;
maskFWHM129 = abs(w_vec - omega129g) < pi * FWHM_max;
[~, Ind129FWHM] = min(abs(yyy .* maskFWHM129 - max(yyy .* maskFWHM129) / 2));
G129g = abs(w_vec(Ind129FWHM) - omega129g);
pos = t < (t(1) + min(max(6 * 2 * pi / omega129g, 0.1), 1 / (3 * G129g)));
R129g = trapz(t(pos), v(pos) .* exp(-1i * omega129g * t(pos))) / (max(t(pos)) - min(t(pos)));
p129g = angle(R129g); 
A129g = 2 * abs(R129g);

maskFWHM131 = abs(w_vec - omega131g) < pi * FWHM_max;
[~, Ind131FWHM] = min(abs(yyy .* maskFWHM131 - max(yyy .* maskFWHM131) / 2));
G131g = abs(w_vec(Ind131FWHM) - omega131g);
pos = t < (t(1) + min(max(6 * 2 * pi / omega131g, 0.1), 1 / (3 * G131g)));
R131g = trapz(t(pos), v(pos) .* exp(-1i * omega131g * t(pos))) / (max(t(pos)) - min(t(pos)));
p131g = angle(R131g); 
A131g = abs(R131g);

if G129g * (t(end) - t(1)) / pi < 1.5
    pos0 = t < (0.1 * t(end) + 0.9 * t(1));
    y0 = trapz(t(pos0), v(pos0) .* exp(-1i * omega129g * t(pos0)));
    pos1 = t > (0.9 * t(end) + 0.1 * t(1));
    y1 = trapz(t(pos1), v(pos1) .* exp(-1i * omega129g * t(pos1)));
    DT = 0.9 * (t(end) - t(1));
    G129g = abs(1 - abs(y1) / abs(y0)) / DT;
end
if G131g * (t(end) - t(1)) / pi < 1.5
    pos0 = t < (0.1 * t(end) + 0.9 * t(1));
    y0 = trapz(t(pos0), v(pos0) .* exp(1i * omega131g * t(pos0)));
    pos1 = t > (0.9 * t(end) + 0.1 * t(1));
    y1 = trapz(t(pos1), v(pos1) .* exp(1i * omega131g * t(pos1)));
    DT = 0.9 * (t(end) - t(1));
    G131g = abs(1 - abs(y1) / abs(y0)) / DT;
end

fo_ = fitoptions('method', 'NonlinearLeastSquares', 'Lower',...
    [-inf, -inf, 0, 0, min(v), 0.9 * omega131g, 0.9 * omega129g, -2*pi, -2*pi],...
    'Upper', [inf, inf, 2 * G131g, 2 * G129g, max(v), 1.1 * omega131g, 1.1 * omega129g, 2 * pi, 2 * pi]);

st_ = [A131g, 2 * A129g, G131g, G129g, mean(v), omega131g, omega129g, p131g, p129g];
set(fo_, 'Startpoint', st_);
ft_ = fittype('a1*exp(-b1*x)*sin(e1*x+f1)+a2*exp(-b2*x)*sin(e2*x+f2)+c',...
    'dependent', {'y'}, 'independent', {'x'}, 'coefficients', {'a1','a2','b1','b2','c','e1','e2','f1','f2'});
try
    if(is_fit_cos == 1)
        ft_ = fittype('a1*exp(-b1*x)*cos(e1*x+f1)+a2*exp(-b2*x)*cos(e2*x+f2)+c', 'dependent', {'y'}, 'independent', {'x'}, 'coefficients', {'a1','a2','b1','b2','c','e1','e2','f1','f2'});
    end
catch
end

x_plot = (t(i_min:i_max) - t(i_min)).'; 
y_plot = v(i_min:i_max).';
[cf_, gof] = fit(x_plot, y_plot, ft_, fo_); % CF = coeffvalues(cf_);

if(to_plot)
    figure(655)
    plot(x_plot, y_plot, 'x-', 'linew', 1.5)
    hold all
    if(length(y_plot) < 1000)
        plot(cf_)
    else
        x_fit = linspace(min(x_plot), max(x_plot), max(length(x_plot), (1e4) + 1));
        y_fit = feval(cf_, x_fit);
        plot(x_fit, y_fit, 'm--', 'linew', 0.5)
    end
    disp(gof)
end
end
