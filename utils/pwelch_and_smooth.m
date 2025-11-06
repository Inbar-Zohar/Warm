function [srPSD,f] = pwelch_and_smooth(input,t_vec,min_t,max_t,use_full)
if ~exist('min_t','var') || isempty(min_t) || isnan(isempty(min_t)); min_t = -Inf; end
if ~exist('max_t','var') || isempty(max_t) || isnan(isempty(max_t)); max_t = Inf; end
if ~exist('use_full','var') || isempty(use_full) || isnan(isempty(use_full)); use_full = 0; end

fs = 1/diff(t_vec(1:2));

inds = t_vec > min_t & t_vec < max_t;
[PSD_full,f_full] = pwelch(input(inds) - mean(input(inds)), [], [], [], fs);
f_full = f_full(2:end-1);
srPSD_full = sqrt(PSD_full(2:end-1));
if use_full
    f = f_full;
    srPSD = srPSD_full;
else
    [f,srPSD] = logarithmic_smooth(f_full, srPSD_full, 'ieee-log');
end
