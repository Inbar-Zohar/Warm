function AD = noise_model_allan_func(tau_s,noise_terms_vals,which_terms)

if ~exist('which_terms','var')
    which_terms = true(1,5);
end
if length(noise_terms_vals) ~= sum(which_terms)
    error('length of noise term vector is %i, but expect total of %i terms', ...
        length(noise_terms_vals), sum(which_terms));
end

tau_hr = tau_s/3600;

powers = -1:0.5:1;
prefactors = sqrt([3 1 2*log(2)/pi 1/3 1/2]);

tau_mat = repmat(tau_hr(:),1,sum(which_terms)) .^ repmat(powers(which_terms),length(tau_hr),1);
AD_terms = repmat(noise_terms_vals(:)'.*prefactors(which_terms),length(tau_hr),1) .* tau_mat;
AD = sqrt(sum(AD_terms.^2,2));