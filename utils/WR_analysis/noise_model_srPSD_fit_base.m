function [noise_term_vals, mdl] = noise_model_srPSD_fit_base(f,srPSD,p0,which_terms)

if ~exist('which_terms','var')
    which_terms = true(1,5);
end
if length(p0) ~= sum(which_terms)
    error('length of initial parameters vector is %i, but expect total of %i terms', ...
        length(p0), sum(which_terms));
end


% 
mdl = fitnlm(f(:),log(srPSD(:)),@(p,f) log(noise_model_srPSD_func(f,abs(p),which_terms)), ...
    p0,'Options',statset('Display','off'));
p_vals = abs(mdl.Coefficients.Estimate)';

noise_term_vals = nan(1,5);
noise_term_vals(which_terms) = p_vals;
