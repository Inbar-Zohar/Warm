function [BWeff,BWeff_err] = get_effective_BW(noise_term_vals,noise_term_errs)

BWeff = (noise_term_vals(1)/noise_term_vals(6))^2;
BWeff_err = BWeff*2*sqrt(sum((noise_term_errs([1 6]) ./ noise_term_vals([1 6])).^2));
