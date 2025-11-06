function full_loop_data = loop_initialize_results_arrays(full_loop_data, params2save)

% Initialize the results struct as nan arrays at the appropriate size for 
% the expected loop.

% if full_loop_data.n_loops > 2 
%     error('No support yet for >2D scans')
if full_loop_data.n_loops == 1
    empty_results_array = nan(1,full_loop_data.n_runs_vec); 
    empty_runs_array = cell(1,full_loop_data.n_runs_vec);
else
    empty_results_array = nan(full_loop_data.n_runs_vec); 
    empty_runs_array = cell(full_loop_data.n_runs_vec);   
end

for ind_param = 1:length(params2save)
    full_loop_data.results.(params2save{ind_param}) = empty_results_array;
end

full_loop_data.individual_runs = empty_runs_array;
