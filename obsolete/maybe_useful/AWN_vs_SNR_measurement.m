%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   The Model:
%                           AWN = a * SNR^-1
%
%   The Experiment Steps:
%       
%       1. Measure the Xe131 signal amplitude from FFT
%       2. Compute the Xe131 noise floor around the peak from its PSD
%       3. Calculate the Xe131 SNR^-1 in units of 1/srHz  
%       4. Measure the Gyroscope PSD
%       5. Extract the AWN from Gyroscope PSD
%       6. Repeat for different SNRs by inserting noise or damaging the
%          signal (detuning the NMR frequency of the 131)
%       7. Fit the model and extract a
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vv = meas.vp_wrld_rot(4, :);
vv_xe_131_noise = meas.vp_wrld_rot(6, :) / 10;
tt = meas.tp_wrld_rot;

[f1, P1, pks, locs, Y] = fft1(tt, vv, 1, 1);

Fs = 1 / mean(diff(tt));
[PSD_dvv, ff] = pwelch(vv - mean(vv), [], [], [], Fs);

ff_size = size(ff);
if ff_size(1) > 1
    PSD_dvv = PSD_dvv';
    ff = ff';
end

low_lim_factor = 2;   % [Hz]
upp_lim_factor = 1.2;   % [Hz]
noise_area_left = all([ff > locs(1) / low_lim_factor; ff < locs(1) / upp_lim_factor]);
noise_area_right = all([ff > locs(1) * upp_lim_factor; ff < locs(1) * low_lim_factor]);
noise_area = any([noise_area_left; noise_area_right]);

linear_fit = fit(ff(noise_area)', sqrt(PSD_dvv(noise_area))', 'poly1');
xe131_noise_eval = feval(linear_fit, locs(1));


% MM = 100;
% NN = length(tt);
% WW = floor((NN - 1) / MM);
% xe131_phase_array = zeros(1, MM);
% for i = 1:MM
%     [f11, ~, ~, ~, YY] = fft1(tt((i-1)*WW+1:i*WW), vv((i-1)*WW+1:i*WW));
%     [~, f1_ind] = min(abs(f11 - locs(1)));
%     f11(f1_ind)
%     xe131_phase_array(i) = angle(YY(f1_ind));
% end
% figure
% plot(1:MM, xe131_phase_array, 'o')
[~, f_ind] = min(abs(f1 - locs(1)));
xe131_phase = angle(Y(f_ind));


figure
loglog(ff, sqrt(PSD_dvv))
hold on
loglog(ff(noise_area), sqrt(PSD_dvv(noise_area)))
loglog(ff(noise_area), feval(linear_fit, ff(noise_area)), 'linew', 3)
loglog(locs(1), xe131_noise_eval, 'o')
loglog(ff(noise_area), mean(vv_xe_131_noise)*ones(size(ff(noise_area))), 'linew', 2)

% get phase from fft
xe131_coherent = pks(1) * cos(2*pi*locs(1) * tt + xe131_phase);

%
vv_no_131 = vv - xe131_coherent;
[PSD_dvv_no_131, ~] = pwelch(vv_no_131 - mean(vv_no_131), [], [], [], Fs);

loglog(ff, sqrt(PSD_dvv_no_131))


figure
plot(tt, vv)
hold on
plot(tt, vv_no_131)



snr = pks(1) / xe131_noise_eval;


[fff, N, T, S, tau, noise_mdl_prm] = get_noise_spectrum_with_fit(meas, prm, prm.estimated_bw);
meas.noise_mdl_prm = noise_mdl_prm;

awn_allan = meas.noise_mdl_prm.allan_PHI_deg_srHz;
awn_psd = meas.noise_mdl_prm.psd_PHI_deg_srHz;

snr * awn_allan