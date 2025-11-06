function [noise_term_vals, noise_term_errs, rmse_vec] = ...
    noise_model_combined_fit(tau,AD,AD_err,f,srPSD,srPSD_err,p0_matrix,...
    rmse_allan,rmse_srPSD,force_keep_terms,allowed_rmse_factor,allowed_contr_factor)

if ~exist('force_keep_terms','var') 
    force_keep_terms = [];
end
if ~exist('allowed_rmse_factor','var') || isempty(allowed_rmse_factor)
    allowed_rmse_factor = 2;
end
if ~exist('allowed_contr_factor','var') || isempty(allowed_contr_factor)
    allowed_contr_factor = 0.25;
end

est_err = (exist('AD_err','var') && ~(isempty(AD_err) || any(isnan(AD_err)))) && ...
    (exist('srPSD','var') && ~(isempty(srPSD) || any(isnan(srPSD))));

inds_matrix = isfinite(p0_matrix(:,1:5));

noise_term_vals = nan(0,6);
noise_term_errs = nan(0,6);
rmse_vec = [];

for ind=1:size(inds_matrix,1)
    [this_noise_term_vals,this_rmse] = noise_model_combined_fit_base(...
        tau,AD,f,srPSD,p0_matrix(ind,[inds_matrix(ind,:) inds_matrix(ind,1)]),...
        rmse_allan,rmse_srPSD,inds_matrix(ind,:));

    force_keep = (~isempty(force_keep_terms) && all(force_keep_terms == inds_matrix(ind,:)));
    this_allan_p0 = this_noise_term_vals(1:5);
    total_AD = noise_model_allan_func(tau,this_allan_p0(inds_matrix(ind,:)),inds_matrix(ind,:));    
    this_srPSD_p0 = this_noise_term_vals([6 2:5]);
    total_srPSD = noise_model_srPSD_func(f,this_srPSD_p0(inds_matrix(ind,:)),inds_matrix(ind,:));

    keep_contr = zeros(1,5);
    for ind_term=find(inds_matrix(ind,:))
        this_AD_term = noise_model_allan_func(...
            tau,this_allan_p0(ind_term),(1:5) == ind_term);
        this_srPSD_term = noise_model_srPSD_func(...
            f,this_srPSD_p0(ind_term),(1:5) == ind_term);
        keep_contr(ind_term) = any(this_AD_term ./ total_AD > allowed_contr_factor) | ...
            any(this_srPSD_term ./ total_srPSD > allowed_contr_factor) ;
    end

    if all(keep_contr(inds_matrix(ind,:))) || force_keep
        noise_term_vals(end+1,:) = this_noise_term_vals;
        rmse_vec(end+1) = this_rmse;

        if est_err           
            this_noise_term_upper_allan = ...
                noise_model_allan_fit_base(tau,AD+AD_err,this_allan_p0(inds_matrix(ind,:)),inds_matrix(ind,:));
            this_noise_term_upper_srPSD = ...
                noise_model_srPSD_fit_base(f,srPSD*(1+srPSD_err(2)),this_srPSD_p0(inds_matrix(ind,:)),inds_matrix(ind,:));

            this_upper_p0 = [this_noise_term_upper_allan(1) ...
                mean([this_noise_term_upper_allan(2:5);this_noise_term_upper_srPSD(2:5)]) ...
                this_noise_term_upper_srPSD(1) ];

            [this_noise_term_upper,~] = noise_model_combined_fit_base(...
                tau,AD+AD_err,f,srPSD*(1+srPSD_err(2)),this_upper_p0(isfinite(this_upper_p0)),...
                rmse_allan,rmse_srPSD,inds_matrix(ind,:));

            noise_term_errs(end+1,:) = abs(this_noise_term_upper - this_noise_term_vals);    
        end
    end
end



[rmse_vec, sorted_inds] = sort(rmse_vec);

noise_term_vals = noise_term_vals(sorted_inds,:);
noise_term_errs = noise_term_errs(sorted_inds,:);

keep_inds = rmse_vec < rmse_vec(1) * allowed_rmse_factor;
if ~isempty(force_keep_terms)
    keep_inds = keep_inds | ismember(isfinite(noise_term_vals(:,1:5)),force_keep_terms,'rows')';
end

noise_term_vals = noise_term_vals(keep_inds,:);
if est_err
    noise_term_errs = noise_term_errs(keep_inds,:);
else
    noise_term_errs = nan(size(noise_term_vals));
end
rmse_vec = rmse_vec(keep_inds);

 
