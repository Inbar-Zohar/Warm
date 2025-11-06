% Extracts into the workspace the relevant variables from the loop data
% struct.

if ~exist('loop_data','var')
    error('Import from loop_data.mat file into variable ''loop_data''');
end

for ind_loop=1:loop_data.n_loops
    this_loop = loop_data.single_loop_data(ind_loop);
    assignin('base',[this_loop.loop_var '_vec'], this_loop.val_vec);
end

res_names = fields(loop_data.results);
for ind_res = 1:length(res_names)
    this_res_name = res_names{ind_res};
    assignin('base',this_res_name, loop_data.results.(this_res_name));
end