function title_str = create_title_str(noise_term_vals, noise_term_errs, type, ...
    prefactors, first_row, tau_range)

if ~exist('noise_term_errs','var') || isempty(noise_term_errs); noise_term_errs = nan*noise_term_vals; end
if ~exist('prefactors','var') || isempty(prefactors); prefactors = [1e-3 1 1 1 1e3]; end
if ~exist('first_row','var'); first_row = ''; end
if ~exist('tau_range','var') || isempty(tau_range); tau_range = [1/100 10*3600]; end

[names, units] = define_variables_for_param_str(type);

which_terms = isfinite(noise_term_vals(1:5));

if strcmp(type,'both')
    title_str_tmp = {'','',''};
    do_both = 1;
    
    [BWeff,BWeff_err] = get_effective_BW(noise_term_vals,noise_term_errs);
    
    [BIeff,BIeff_err,is_bound,~] = get_effective_BI(noise_term_vals,noise_term_errs,tau_range);
    
    which_terms = logical([which_terms([1 1 1 2 3]) 1 which_terms([4 5])]);
    noise_term_vals = [noise_term_vals([1 6]) BWeff ...
        noise_term_vals([2 3]) BIeff noise_term_vals([4 5])];
    noise_term_errs = [noise_term_errs([1 6]) BWeff_err ...
        noise_term_errs([2 3]) BIeff_err noise_term_errs([4 5])];
    prefactors = [prefactors([1 1]) 1 prefactors([2 3 3 4 5])];
else
    title_str_tmp = {'',''};
    do_both = 0;
end



for ind = find(which_terms)
    this_str = create_variable_str(names{ind}, ...
        noise_term_vals(ind), noise_term_errs(ind), prefactors(ind), units{ind});
    if do_both && ind == 6 && is_bound
        this_str = strrep(this_str,'=','\leq');
    end
    if do_both
        this_row = floor((ind-1)/3) + 1;
    else
        this_row = 1 + (ind>2) ;
    end
    title_str_tmp{this_row} = [title_str_tmp{this_row} '; ', this_str];
end

title_str = {};
for row=1:length(title_str_tmp)
    if ~isempty(title_str_tmp{row})
        title_str{end+1} = title_str_tmp{row}(3:end);
    end
end

if ~isempty(first_row)
    title_str = [first_row,title_str];
end