function [noise_term_vals,rmse] = noise_model_combined_fit_base(tau,AD,f,srPSD,p0,...
    rmse_allan,rmse_srPSD,which_terms)

if ~exist('which_terms','var')
    which_terms = true(1,5);
end
if length(p0) ~= sum(which_terms) + which_terms(1)
    error('length of initial parameters vector is %i, but expect total of %i terms', ...
        length(p0), sum(which_terms) + which_terms(1));
end

n_terms = sum(which_terms);

allan_func = @(p) noise_model_allan_func(tau,p(1:n_terms),which_terms);
if which_terms(1)
    srPSD_func = @(p) noise_model_srPSD_func(f,p([n_terms+1 2:n_terms]),which_terms);
else
    srPSD_func = @(p) noise_model_srPSD_func(f,p(1:n_terms),which_terms);
end
fmin_func = @(p) (1/2 * ( mean( log(AD(:)./allan_func(p)).^2 )/rmse_allan.^2 + ...
    mean( log(srPSD(:)./srPSD_func(p)).^2 )/rmse_srPSD.^2 ));

options = optimoptions('fminunc');
options.MaxFunctionEvaluations = 3e3;
options.Display = 'none';

[p_vals,rmse,~,~,~,~] = fminunc(fmin_func, p0,options);

which_terms_comb = [which_terms which_terms(1)];

noise_term_vals = nan(1,6);
noise_term_vals(which_terms_comb) = abs(p_vals);
