function legend_str = create_legend_str(noise_term_vals, noise_term_errs, type, prefactors)

[~, units] = define_variables_for_param_str(type, prefactors);

which_terms = isfinite(noise_term_vals);

legend_str = '';

for ind = find(which_terms)
    this_str = create_variable_str('', ...
        noise_term_vals(ind), noise_term_errs(ind), prefactors(ind), units{ind});
    
%     prefactor = prefactors(ind);
%     if isfinite(noise_term_errs(ind))
%         this_str = sprintf('%s %s', ...
%             numerr2str(noise_term_vals(ind)/prefactor, noise_term_errs(ind)/prefactor), units{ind});
%     else
%         this_str = sprintf('%s %s', ...
%             num2str(noise_term_vals(ind)/prefactor,2), units{ind});
%     end
    
    legend_str = [legend_str '; ', this_str];
end

if ~isempty(legend_str)
    legend_str = legend_str(3:end);
end