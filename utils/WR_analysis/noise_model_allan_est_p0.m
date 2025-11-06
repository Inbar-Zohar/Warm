function p0 = noise_model_allan_est_p0(tau_s,AD)

powers = -1:0.5:1;
prefactors = sqrt([3 1 1 1/3 1/2]);
tau_hr = tau_s/3600;

tau_mat = repmat(tau_hr(:),1,length(powers)) .^ repmat(powers,length(tau_hr),1);
X = repmat(prefactors,length(tau_hr),1) .* tau_mat;

warning('off','MATLAB:nearlySingularMatrix')
p0 = abs((X'*X + 5e-3*eye(5))\X'*(AD(:)));
warning('on','MATLAB:nearlySingularMatrix')
