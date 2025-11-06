function [p0_matrix,bw_est_from_fits] = ...
    create_p0_matrix_for_combined_fit(allan_noise_term_vals, srPSD_noise_term_vals, p00)

if ~exist('p00','var') || isempty(p00)
    p00 = zeros(1,5);
end

bw_est_from_fits = (nanmean(allan_noise_term_vals(:,1)) ./  ...
    nanmean(srPSD_noise_term_vals(:,1)))^2;
if isnan(bw_est_from_fits); bw_est_from_fits = 1; end;

if 0
    terms_matrix = unique(isfinite([allan_noise_term_vals ; srPSD_noise_term_vals]),'rows');
else
    [rows_allan, ~] = size(allan_noise_term_vals);
    [rows_srPSD, ~] = size(srPSD_noise_term_vals);
    combinedOrRows = reshape(isfinite(allan_noise_term_vals), [rows_allan, 1, 5]) ...
        | reshape(isfinite(srPSD_noise_term_vals), [1, rows_srPSD, 5]);
    terms_matrix = unique([isfinite(allan_noise_term_vals);
        isfinite(srPSD_noise_term_vals);
        reshape(combinedOrRows, [rows_allan * rows_srPSD, 5])],'rows');
end
p0_matrix = nan(size(terms_matrix,1),6);
for ind=1:size(terms_matrix,1)
    this_terms = terms_matrix(ind,:);
    this_allan_ind = find(ismember(isfinite(allan_noise_term_vals),this_terms,'rows'));
    this_psd_ind = find(ismember(isfinite(srPSD_noise_term_vals),this_terms,'rows'));
    if ~isempty(this_allan_ind) && ~isempty(this_psd_ind)
        p0_matrix(ind,1) = allan_noise_term_vals(this_allan_ind,1);
        p0_matrix(ind,6) = srPSD_noise_term_vals(this_psd_ind,1);
        p0_matrix(ind,2:5) = mean([allan_noise_term_vals(this_allan_ind,2:5);...
            srPSD_noise_term_vals(this_psd_ind,2:5)]);
    elseif ~isempty(this_allan_ind) && isempty(this_psd_ind)
        p0_matrix(ind,1:5) = allan_noise_term_vals(this_allan_ind,:);
        p0_matrix(ind,6) = p0_matrix(ind,1) / sqrt(bw_est_from_fits);
    elseif isempty(this_allan_ind) && ~isempty(this_psd_ind)
        p0_matrix(ind,[6 2:5]) = srPSD_noise_term_vals(this_psd_ind,:);
        p0_matrix(ind,1) = p0_matrix(ind,6) * sqrt(bw_est_from_fits);
    else
        all_awn = nanmean([allan_noise_term_vals(:,1);srPSD_noise_term_vals(:,1)*sqrt(bw_est_from_fits)]);

        this_p0 = [all_awn,...
            nanmean([allan_noise_term_vals(:,2:5);srPSD_noise_term_vals(:,2:5)]),...
            all_awn/sqrt(bw_est_from_fits)];
        this_p0(~this_terms) = nan;
        if any(~isfinite(this_p0(this_terms([1:5 1]))))
            keyboard
            this_p0 = [allan_noise_term_vals(:,[1:5 1]); srPSD_noise_term_vals(:,[1:5 1])];
            error('Fucked up somehow in the p0 defintions :(')
        end
        p0_matrix(ind,:) = this_p0;
    end
end
return
%%
% inds_matrix_all = create_all_combinations;
% inds_matrix_allan = isfinite(allan_noise_term_vals);
% inds_matrix_srPSD = isfinite(srPSD_noise_term_vals);
% 
% for ind=1:size(inds_matrix_all,1)
%     this_terms = inds_matrix_all(ind,:);
%     this_terms = terms_matrix(ind,:);
%     this_allan_ind = find(ismember(isfinite(allan_noise_term_vals),this_terms,'rows'));
%     this_psd_ind = find(ismember(isfinite(srPSD_noise_term_vals),this_terms,'rows'));
%     this_p0_vec = nan(1,6);
%     this_p0_vec
%     if ~isempty(this_allan_ind) && ~isempty(this_psd_ind)
%         p0_matrix(ind,1) = allan_noise_term_vals(this_allan_ind,1);
%         p0_matrix(ind,6) = srPSD_noise_term_vals(this_psd_ind,1);
%         p0_matrix(ind,2:5) = mean([allan_noise_term_vals(this_allan_ind,2:5);...
%             srPSD_noise_term_vals(this_psd_ind,2:5)]);
%     elseif ~isempty(this_allan_ind) && isempty(this_psd_ind)
%         p0_matrix(ind,1:5) = allan_noise_term_vals(this_allan_ind,:);
%         p0_matrix(ind,6) = p0_matrix(ind,1) / sqrt(bw_est_from_fits);
%     elseif  isempty(this_allan_ind) && ~isempty(this_psd_ind)
%         p0_matrix(ind,[6 2:5]) = srPSD_noise_term_vals(this_psd_ind,:);
%         p0_matrix(ind,1) = p0_matrix(ind,6) * sqrt(bw_est_from_fits);
%     else
%         matching_rows_allan = [];
%         matching_rows_srPSD = [];
%         
%         for i = 1:size(inds_matrix_allan, 1)
%             if all(inds_matrix_allan(i, inds_matrix_allan(i, :)) <= this_terms(inds_matrix_allan(i, :)))
%                 matching_rows_allan = [matching_rows_allan ind];
%             end
%         end
%         for i = 1:size(inds_matrix_allan, 1)
%             if all(inds_matrix_allan(i, inds_matrix_allan(i, :)) <= this_terms(inds_matrix_allan(i, :)))
%                 matching_rows_allan = [matching_rows_allan ind];
%             end
%         end
% 
% %         if ~(isempty(matching_rows_allan) && matching_rows_allan
%         if ~isempty(matching_rows)
%             this_terms_full = this_terms([1:5 1]);
%             this_p0_vec = nan(1,6);
%             this_p0_vec(this_terms_full) == p00(this_terms_full);
%         end
% 
%         
%     end
% end
% return
