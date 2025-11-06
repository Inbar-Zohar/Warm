%%ADEV

function [sigy, tau] = ADEV(Sig, fs)

dT = 1 / fs; %mean sampling time [sec]

N = length(Sig);
%% outliers removal


%% ADEV
%tau=[1 5 10 50 100 500  1e3 5e3 1e4 5e4]*dT; %[sec]
m = 300;
tauP = dT * (logspace(log10(1), log10(N-1), m));
tau = floor(tauP ./ tauP(1)) * dT;

for k = 1:length(tau)
    k;
    Sig2 = Sig;
    j = tau(k)/dT;
    
    % binning
    res = mod(N, j);
    Sig2(end-res+1:end) = [];
    Sig3 = reshape(Sig2, int32(j), int32(length(Sig2) / j));
    Sig3 = sum(Sig3) / j;
    
    ds3 = diff(Sig3);
    N3 = 2 * (length(Sig3) - 1);%normalization factor per tau
    AVAR(k) = 1 / N3 * sum(ds3.^2);
    ADV(k) = sqrt(AVAR(k));
    
    clear  Sig2 Sig3 N3 ds3 res
end

sigy = ADV;

end