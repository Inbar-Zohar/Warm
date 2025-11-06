function P = NMRG_phase_freq_transfer_func(s, T2, varargin)

% returns the phase-frequency transfer function on an NMRG according to
% PRAppl 5, 014042 (2024) in units of 1/s

p = inputParser;

ispositive = @(x) isnumeric(x) && x>0;
addRequired(p,'s',@isnumeric);
addRequired(p,'T2',@(x)ispositive(x));

addParameter(p,'kP',10*max(abs(s)),@isnumeric);
addParameter(p,'kI',10*max(abs(s))/T2,@isnumeric);

isfilter = @(x) isstruct(x) && isfield(x,'func') && isfield(x,'params') && ...
    isa(x.func,'function_handle') && isnumeric(x.params);
default_filter = struct('func',@(x,params) 1, 'params', []);
addParameter(p,'filter_num',default_filter,@(x)isfilter(x));
addParameter(p,'filter_denom',default_filter,@(x)isfilter(x));

parse(p, s, T2, varargin{:});
struct2var(p.Results);

C = kP + kI./s;
G = 1./ (1./T2 + s); 
F_num = filter_num.func(s,filter_num.params);
F_denom = filter_denom.func(s,filter_denom.params);
% Total transfer functions
P = C .* F_num ./ (1 + C .* G .* F_denom);


