function srPSD = noise_model_srPSD_func(f,noise_terms_vals,which_terms)

if ~exist('which_terms','var')
    which_terms = true(1,5);
end
if length(noise_terms_vals) ~= sum(which_terms)
    error('length of noise term vector is %i, but expect total of %i terms', ...
        length(noise_terms_vals), sum(which_terms));
end

omega = 2*pi*f * 3600;
powers = [1 0 -0.5 -1 -1.5];
prefactors = sqrt(2)*[1/60 1 1 1 1];

omega_mat = repmat(omega(:),1,sum(which_terms)) .^ repmat(powers(which_terms),length(omega),1);
PSD_terms = repmat(noise_terms_vals(:)'.*prefactors(which_terms),length(omega),1) .* omega_mat;
srPSD = sqrt(sum(PSD_terms.^2,2));
