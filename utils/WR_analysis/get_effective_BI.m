function [BIeff,BIeff_err,is_bound,tau_min] = get_effective_BI(noise_term_vals,noise_term_errs,tau_range)

noise_term_vals = noise_term_vals(1:5);
noise_term_errs = noise_term_errs(1:5);
which_terms = isfinite(noise_term_vals);
if find(which_terms, 1, 'last' ) == 3
    BIeff = noise_term_vals(3);
    is_bound = 0;
    BIeff_err = noise_term_errs(3);

    prefactors = sqrt([3 1 2*log(2)/pi]);
    allan_terms = prefactors .* noise_term_vals(1:3);
    tau_min = nanmax([allan_terms(1)/allan_terms(3),...
        (allan_terms(2)/allan_terms(3))^2]*3600);
    
else
    [tau_min,AD_min] = fminbnd(...
        @(tau) noise_model_allan_func(tau,noise_term_vals(which_terms),which_terms),...
        tau_range(1),tau_range(2));
    if tau_min > 0.95*tau_range(2) || tau_min < 0.05*tau_range(1)
        is_bound = 1;
    else
        is_bound = 0;
    end
    BI_factor = sqrt(2*log(2)/pi);
    BIeff = AD_min / BI_factor;
    noise_term_vals_upper = noise_term_vals + noise_term_errs;
    BIeff_err = (...
        noise_model_allan_func(tau_min,noise_term_vals_upper(which_terms),which_terms) - AD_min) ...
        / BI_factor;
end
