function full_loop_data = loop_update_loop_results(full_loop_data,meas,folder_path)

% Extracts the desired results from the meas struct to the full_loop_data
% results struct.
% On the first run, check if the desired parameters to save have been
% defined. If not, add all scalar fields from this meas struct.

params2save = fields(full_loop_data.results);

if isempty(params2save) && full_loop_data.ind_run == 1
    all_meas_params = fields(meas);
    inds2keep = zeros(size(all_meas_params));
    for ind_param = 1:length(all_meas_params)
        this_param = all_meas_params{ind_param};
        inds2keep(ind_param) = ~strcmp('DateTime_var',this_param) & ...
            ~strcmp('DateNumber',this_param) & ...
            numel(meas.(this_param)) == 1;
    end
    params2save = all_meas_params(logical(inds2keep));
    
    full_loop_data = loop_initialize_results_arrays(full_loop_data,params2save);
end

if full_loop_data.ind_run == 1
    full_loop_data.DateTime_var = meas.DateTime_var;
    full_loop_data.DateNumber = meas.DateNumber;
end

for ind_param = 1:length(params2save)
    this_param = params2save{ind_param};
    if isfield(meas,this_param) && numel(meas.(this_param)) == 1
        full_loop_data.results.(this_param)(full_loop_data.ind_run) = ...
            meas.(this_param);
    end
end

full_loop_data.individual_runs{full_loop_data.ind_run} = folder_path;