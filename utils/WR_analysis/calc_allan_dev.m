% copied from http://www.nxp.com/assets/documents/data/en/application-notes/AN5087.pdf
function [T, sigma, sigma_err] = calc_allan_dev(omega, fs, pts)

omega = omega(:);
[N,M] = size(omega); % figure out how big the output data set is
n = 2.^(floor(log2(N/2)))'; % determine largest bin size
m = unique(ceil(logspace(0,log10(n),pts)))'; % create log spaced vector average factor
% 	m = m(1:end-1);
t0 = 1/fs; % t0 = sample interval
T = m*t0; % T = length of time for each cluster
theta = cumsum(omega)/fs; % integration of samples over time to obtain output angle ?
sigma2 = zeros(length(T),M); % array of dimensions (cluster periods) X (#variables)
for i=1:length(m) % loop over the various cluster sizes
    for k=1:N-2*m(i) % implements the summation in the AV equation
        sigma2(i,:) = sigma2(i,:) + (theta(k+2*m(i),:) - 2*theta(k+m(i),:) + theta(k,:)).^2;
    end
end
sigma2 = sigma2./repmat((2*T.^2.*(N-2*m)),1,M);
sigma = sqrt(sigma2);
sigma_err = sigma./sqrt(n./m); % see page 37 of https://nvlpubs.nist.gov/nistpubs/Legacy/SP/nistspecialpublication1065.pdf


