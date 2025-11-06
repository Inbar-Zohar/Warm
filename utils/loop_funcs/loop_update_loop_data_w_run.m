function full_loop_data = loop_update_loop_data_w_run(full_loop_data, ind_run)

% Add to the full_loop_data the information about the current run. 

full_loop_data.ind_run = ind_run;
if full_loop_data.n_loops == 1
    full_loop_data.sub_run = ind_run;
else
%     [ind1, ind2] = ind2sub(full_loop_data.n_runs_vec,ind_run);
%     full_loop_data.sub_run = [ind1, ind2];
    out = cell(size(full_loop_data.n_runs_vec));
    [out{:}] = ind2sub(full_loop_data.n_runs_vec,ind_run);
    full_loop_data.sub_run = [out{:}];
end

full_loop_data.val_run = nan(1,full_loop_data.n_loops);
for ind_single = 1:length(full_loop_data.sub_run)
    this_loop = full_loop_data.single_loop_data(ind_single);
    full_loop_data.val_run(ind_single) = this_loop.val_vec(...
        full_loop_data.sub_run(ind_single));
end