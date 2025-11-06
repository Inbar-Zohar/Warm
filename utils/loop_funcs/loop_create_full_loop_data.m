function full_loop_data = loop_create_full_loop_data(single_loop_data, params2save)

% Create the full_loop_data struct from the user-defined single_loop_data
% and params2save.

for ind_loop = 1:length(single_loop_data)
    this_lp = single_loop_data(ind_loop);
    if isfield(this_lp,'unit') && ~isempty(this_lp.unit)
        loop_var_display =  [this_lp.loop_var ' (' this_lp.unit ')'];
    else
        this_lp.unit = '';
        loop_var_display = this_lp.loop_var;
    end
    single_loop_data(ind_loop).loop_var_display = loop_var_display;

    if ~isfield(this_lp,'calib')
        single_loop_data(ind_loop).calib = 1;
    end
    
    single_loop_data(ind_loop).n_runs = ...
        length(this_lp.val_vec);
end

full_loop_data.single_loop_data = single_loop_data;
full_loop_data.results = struct;
full_loop_data.n_runs_vec = [single_loop_data.n_runs];
full_loop_data.n_loops = length(full_loop_data.n_runs_vec);
full_loop_data.n_runs_total = prod(full_loop_data.n_runs_vec);

full_loop_data = loop_initialize_results_arrays(full_loop_data,params2save);

