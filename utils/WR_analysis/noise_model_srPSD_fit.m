function [noise_term_vals, noise_term_errs, mdls, rmse] = ...
    noise_model_srPSD_fit(f,srPSD,srPSD_err,p0,force_keep_terms,allowed_rmse_factor,allowed_contr_factor)

if ~exist('force_keep_terms','var') 
    force_keep_terms = [];
end
if ~exist('allowed_rmse_factor','var') || isempty(allowed_rmse_factor)
    allowed_rmse_factor = 2;
end
if ~exist('allowed_contr_factor','var') || isempty(allowed_contr_factor)
    allowed_contr_factor = 0.25;
end

noise_term_vals = nan(0,5);
noise_term_errs = nan(0,5);
mdls = {};
rmse = [];

%% create all relevant subsets of parameters (at least one decaying and one stable/increasing)
array1 = [1 0; 0 1; 1 1];
array2 = [0 0 0; 1 0 0; 0 1 0; 0 0 1; 1 1 0; 1 0 1; 0 1 1; 1 1 1];
repeated_array1 = repmat(array1, size(array2, 1), 1);
repeated_array2 = kron(array2, ones(size(array1, 1), 1));
combinations = [repeated_array1, repeated_array2];
inds_matrix = logical(combinations);

%% do fits on all subsets
for ind=1:size(inds_matrix,1)
    [this_noise_term_vals, this_mdl] = noise_model_srPSD_fit_base(...
        f,srPSD,p0(inds_matrix(ind,:)),inds_matrix(ind,:));
    force_keep = (~isempty(force_keep_terms) && all(force_keep_terms == inds_matrix(ind,:)));

    if all(this_mdl.Coefficients.pValue < 0.05) || force_keep
        total_srPSD = noise_model_srPSD_func(...
            f,this_noise_term_vals(inds_matrix(ind,:)),inds_matrix(ind,:));
        keep_contr = zeros(1,5);        
        for ind_term=find(inds_matrix(ind,:))
            this_srPSD_term = noise_model_srPSD_func(...
                f,this_noise_term_vals(ind_term),(1:5) == ind_term);
            keep_contr(ind_term) = ...
                any(this_srPSD_term ./ total_srPSD > allowed_contr_factor);
        end
        if all(keep_contr(inds_matrix(ind,:))) || force_keep
            noise_term_vals(end+1,:) = this_noise_term_vals;
            mdls{end+1} = this_mdl;
            rmse(end+1) = this_mdl.RMSE;
        end
    end
end

keep_inds = (rmse < allowed_rmse_factor*min(rmse));
if ~isempty(force_keep_terms)
    keep_inds = keep_inds | ismember(isfinite(noise_term_vals),force_keep_terms,'rows')';
end
noise_term_vals = noise_term_vals(keep_inds,:);
mdls = mdls(keep_inds);
rmse = rmse(keep_inds);

[rmse, sorted_inds] = sort(rmse);
noise_term_vals = noise_term_vals(sorted_inds,:);
mdls = mdls(sorted_inds);

if exist('srPSD','var') && ~(isempty(srPSD) || any(isnan(srPSD)))
    for ind=1:size(noise_term_vals,1)
        which_terms = isfinite(noise_term_vals(ind,:));
        this_noise_term_w_errs = noise_model_srPSD_fit_base(...
            f,srPSD*(1+srPSD_err(2)),noise_term_vals(ind,which_terms),which_terms);
        noise_term_errs(ind,:) = abs(this_noise_term_w_errs - noise_term_vals(ind,:));
    end
else
    noise_term_errs = nan(size(noise_term_vals));    
end
