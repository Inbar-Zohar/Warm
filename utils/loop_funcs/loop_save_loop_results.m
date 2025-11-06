function full_loop_data = loop_save_loop_results(full_loop_data,new_folder,exp_folder_path)

% Saves full_loop_data both in the current experiment file, and in the
% dedicated folder for the loop (where it'll be overwritten by the later
% iterations). If the dedicated folder does not exist, it's created.


if new_folder || ~isfield(full_loop_data, 'loop_folder_path') || ...
        ~isfolder(full_loop_data.loop_folder_path)
    loop_folder_path = exp_folder_path;
    if loop_folder_path(end) == '/' || loop_folder_path(end) == '\'
        loop_folder_path = loop_folder_path(1:end-1);
    end
    loop_folder_path = [loop_folder_path '_loop'];
    full_loop_data.loop_folder_path = loop_folder_path;

    mkdir(loop_folder_path)
    copyfile('utils\loop_funcs\default_loop_analysis.m',fullfile(loop_folder_path,'analysis.m'))
end

save(fullfile(full_loop_data.loop_folder_path, 'loop_data'), '-struct', 'full_loop_data');
save(fullfile(exp_folder_path, 'loop_data'), '-struct', 'full_loop_data');
