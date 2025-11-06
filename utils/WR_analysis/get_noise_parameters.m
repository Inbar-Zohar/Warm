function [tau_range, f_range,...
    noise_term_vals, noise_term_errs, ...
    allan_noise_term_vals, allan_noise_term_errs, ...
    srPSD_noise_term_vals, srPSD_noise_term_errs, ...
    rmse_vec, allan_rmse, srPSD_rmse] = ...
    get_noise_parameters(tau, AD, AD_err, f, srPSD, srPSD_err, ...
    bw_est, force_keep_terms)
% Performs analysis of the Allan deviation and (sqrt of) power spectral 
% density to find all sets of noise parameters consistent with the either 
% one or both curves. 
%
% Input parameters:
%   tau, AD, AD_err: Allan deviation data (integration time, ADEV, ADEV
%       uncertainty). tau in seconds, AD and AD_err in deg/hr.
%   f, srPSD, srPSD_err: srPSD data (frequency, srPSD, srPSD uncertainty).
%       f in Hz, srPSD in deg/srhr, srPSD_err a two-element vector
%       of relative bottom and top uncertainty.
%   be_est: the estimated bandwidth of the measurement. Used to determine
%       the fitting range.
%   force_keep_terms: optional 1x5 boolean vector defining a preferred set
%       of noise terms that should be kept.
%
% Output variables:
%   tau_range, f_range: the ranges over which the model was fitted.
%   noise_term_vals, noise_term_errs: the results of the combined fit to
%       both noise curve.
%   allan_noise_term_vals, allan_noise_term_errs: the results of the fit
%       only to the Allan.
%   srPSD_noise_term_vals, srPSD_noise_term_errs: the results of the fit
%       only to the srPSD.
%   rmse_vec, allan_rmse, srPSD_rmse: the RMSE of the sets of fits
%       (combined, only Allan, only srPSD).

%% Find the ranges of the curves to carry out the fits
[allan_min_ind,allan_max_ind,srPSD_min_ind,srPSD_max_ind] = ...
    get_fit_range(tau,f,bw_est);
allan_inds = allan_min_ind:allan_max_ind;
srPSD_inds = srPSD_min_ind:srPSD_max_ind;
tau_range = tau([allan_min_ind allan_max_ind]);
f_range = f([srPSD_min_ind srPSD_max_ind]);

%% Find noise models consistent with the Allan
[allan_noise_term_vals, allan_noise_term_errs, ~, allan_rmse] = ...
    noise_model_allan_fit(tau(allan_inds),AD(allan_inds),AD_err(allan_inds),force_keep_terms);

%% Find noise models consistent with the PSD
p0 = noise_model_allan_est_p0(tau(allan_inds),AD(allan_inds));
p0(isfinite(allan_noise_term_vals(1,:))) = allan_noise_term_vals(1,isfinite(allan_noise_term_vals(1,:)));
[srPSD_noise_term_vals, srPSD_noise_term_errs, ~, srPSD_rmse] = ...
    noise_model_srPSD_fit(f(srPSD_inds),srPSD(srPSD_inds),srPSD_err,p0,force_keep_terms);

%% Find noise models consistent with both
[p0_matrix,~] = ...
    create_p0_matrix_for_combined_fit(allan_noise_term_vals, srPSD_noise_term_vals);
[noise_term_vals, noise_term_errs, rmse_vec] = ...
    noise_model_combined_fit(tau(allan_inds),AD(allan_inds),AD_err(allan_inds),...
    f(srPSD_inds),srPSD(srPSD_inds),srPSD_err,p0_matrix,...
    allan_rmse(1),srPSD_rmse(1),force_keep_terms,3);
