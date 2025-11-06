function [tau, AD, AD_err, f, srPSD, srPSD_err, f_full, srPSD_full] = ...
    get_allan_and_psd_from_data(Omega_deg_hr, fs, allan_opts, psd_opts)

% Performs analysis of a time trace of a rotation signal returning the the 
% Allan deviation and (sqrt of) power spectral density. 
% Input parameters: time trace of the rotation signal in deg/hr, sampling 
% frequnecy in Hz, and optional strucs with analysis options. 
% For the allan deviation, the calculation is straight-forward based on the
% the internal calc_allan_dev function. The one option is determining
% number of points. 
% For the PSD, following the calculation using built-in function pwelch,
% there's an optional (on by default) stage of downsampling and filtering,
% to get a smooth curve sampled logarithmically. Options include a flag
% determining whether to execute that downsampling and smoothing, and
% parameters related to it. First and last points out of pwelch have higher
% uncertainty and are dropped.
%
% Input parameters:
%   Omega_deg_hr: rotation signal to be analyzed, in deg/hr. 
%   fs: sampling frequency in 1/s. 
%   allan_opts: optional struct with Allan analysis options. Specifically,
%       can determine the number of points requested on the Allan curve.
%   psd_opts: optional struct with PSD analysis options. Specifically, it
%       contains the flag and options for the logarithmic downsampling. 

% Output variables:
%   tau: Allan deviation integration time, in s
%   AD: Allan deviation values, in deg/hr
%   AD_err: Allan deviation base values (additive, same length as AD), in deg/hr
%   f: PSD frequency after downsampling and smoothing, in Hz 
%   srPSD: sqrt(PSD) after downsampling and smoothing, in deg/sqrt(hr)
%   srPSD_err: uncertainty of sqrt(PSD) in 67% CI, first/second term is
%       multiplicative factor for negative/positive error
%   f_full: PSD frequency without downsampling and smoothing, in Hz 
%   srPSD_full: sqrt(PSD) without downsampling and smoothing, in deg/sqrt(hr)


%% Allan
if ~exist('allan_opts','var') || isempty(allan_opts); allan_opts = struct; end
if ~isfield(allan_opts,'Npoints'); allan_opts.Npoints = 100; end

[tau, AD, AD_err] = calc_allan_dev(Omega_deg_hr, fs, allan_opts.Npoints);

% Last point uncertainty is equal to last point, decreasing a bit to allow
% display on a log scale
AD_err(end) = 0.99*AD_err(end); 

%% PSD
if ~exist('psd_opts','var') || isempty(psd_opts); psd_opts = struct; end
if ~isfield(psd_opts,'downsample_and_filter'); psd_opts.downsample_and_filter = 1; end
if psd_opts.downsample_and_filter
    if ~isfield(psd_opts,'method'); psd_opts.method = 'ieee-log'; end
    if ~isfield(psd_opts,'smooth_movmean'); psd_opts.smooth_movmean = 20; end
    if ~isfield(psd_opts,'output_length'); psd_opts.output_length = 100; end
    if ~isfield(psd_opts,'interp_method'); psd_opts.interp_method = 'pchip'; end
end

[dOmega_PSD_deg2_hr2_hz,f,dOmega_PSD_deg2_hr2_hz_confint] = ...
    pwelch(Omega_deg_hr - mean(Omega_deg_hr), [], [], [], fs, 'ConfidenceLevel',0.67);   
f = f(2:end-1);
srPSD = sqrt(dOmega_PSD_deg2_hr2_hz(2:end-1) / 3600);

f_full = f;
srPSD_full = srPSD;
if psd_opts.downsample_and_filter
    [f,srPSD] = logarithmic_smooth(f_full, srPSD_full, psd_opts.method, psd_opts);
end

srPSD_err = abs(mean(sqrt( ...
    dOmega_PSD_deg2_hr2_hz_confint(2:end-1,:) ./ dOmega_PSD_deg2_hr2_hz(2:end-1)  ...
    )) - 1);