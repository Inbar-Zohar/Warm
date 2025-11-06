function loop_description = loop_parse_loop_data(loop_data)

% Create loop description string from loop data struct, e.g. for the
% logger.

loop_description = 'Run is part of a loop.';
try
    nruns_str = [num2str(loop_data.ind_run) '/' num2str(loop_data.n_runs_total)];

    var_str = '';
    for ind_loop = 1:loop_data.n_loops
        this_var_str = [loop_data.single_loop_data(ind_loop).loop_var ...
            ' -> ' num2str(loop_data.val_run(ind_loop))];
        if ~isempty(loop_data.single_loop_data(ind_loop).unit)
            this_var_str = [this_var_str ' ' loop_data.single_loop_data(ind_loop).unit];
        end
        var_str = [var_str this_var_str ', '];
    end
    var_str = var_str(1:end-2);
    
    loop_description = [loop_description(1:end-1) ' : ' nruns_str ' ; ' var_str];
catch
%     loop_description = [loop_description '. Failed parsing loop data.'];   
end