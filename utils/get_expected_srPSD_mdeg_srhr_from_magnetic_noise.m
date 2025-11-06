function srPSD_vec_mdeg_sqrthr = get_expected_srPSD_mdeg_srhr_from_magnetic_noise(f_vec,...
    Gamma129, Gamma131, kP, kI129, kI131, k_f, tau_f, noise_floor_fTsqrtHz, EMF129_uG, EMF131_uG)


% Filter functions
Gfunc = @(s,Gamma) 1./(s + Gamma);
Cfunf = @(s,kP,kI) kP + kI./s;
Ffunc = @(s,k_f,tau_f) 1./(1+tau_f*s).^k_f;
Hfunc = @(s,Gamma,kP,kI,k_f,tau_f) ...
    Gfunc(s,Gamma).*Cfunf(s,kP,kI).*Ffunc(s,k_f,tau_f) ./ ...
    (1 + Gfunc(s,Gamma).*Cfunf(s,kP,kI).*Ffunc(s,k_f,tau_f));
Pfunc = @(s,Gamma,kP,kI,k_f,tau_f) ...
    Cfunf(s,kP,kI).*Ffunc(s,k_f,tau_f) ./ ...
    (1 + Gfunc(s,Gamma).*Cfunf(s,kP,kI).*Ffunc(s,k_f,tau_f));


R = 3.3734;


% Magnetic noise floor converted to mdeg
S_phi_w_129_mdeg2_Hz = (180e3 / pi * noise_floor_fTsqrtHz / (EMF129_uG*1e-6*1e-4*1e15))^2;
S_phi_w_131_mdeg2_Hz = (180e3 / pi * noise_floor_fTsqrtHz / (EMF131_uG*1e-6*1e-4*1e15))^2;

% 
% PSD_vec_mdeg2_Hz_129 = ...
%     abs(Pfunc(1i*2*pi*f_vec,Gamma129,kP,kI129,k_f,tau_f)).^2/(1+R)^2 * (S_phi_w_129_mdeg2_Hz);
% PSD_vec_mdeg2_Hz_131 = ...
%     abs(Pfunc(1i*2*pi*f_vec,Gamma131,kP,kI131,k_f,tau_f)).^2/(1+R)^2*R^2 * (S_phi_w_131_mdeg2_Hz);

srPSD_vec_mdeg_sqrthr = sqrt( 3600 * (...
    abs(Pfunc(1i*2*pi*f_vec,Gamma129,kP,kI129,k_f,tau_f)).^2/(1+R)^2 * S_phi_w_129_mdeg2_Hz + ...
    abs(Pfunc(1i*2*pi*f_vec,Gamma131,kP,kI131,k_f,tau_f)).^2/(1+R)^2*R^2 * S_phi_w_131_mdeg2_Hz) );