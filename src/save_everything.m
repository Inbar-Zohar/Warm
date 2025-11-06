% ========================================================================
% Save raw data, computed variables and all open plots
% This script saves the main script, variables and open figures
% to the backup experiment folder
% ========================================================================

% create new folder for experiment
set_folder_for_save


% save all figures created in the code
if ~isfield(meas,'figs')
    fig_handles = [];
else
    fig_handles = meas.figs;
    meas = rmfield(meas,'figs');
end

if flg.save_all_open_figs
    fig_handles = findall(0, 'Type', 'figure');
end


% save raw data and variables
save(fullfile(folder_path, 'experiment_raw_data'), 'meas', 'prm', 'flg');

if not(isempty(fig_handles))
    for i = 1:length(fig_handles)
        figname = num2str(fig_handles(i).Number);
        try
            savefig(fig_handles(i), fullfile(folder_path, 'plots', figname))
            print(fig_handles(i), fullfile(folder_path, 'plots', figname), '-dpng','-r200');
        catch e
            msg = ['A problem occured while saving figure ' num2str(figname) '. Other data is saved. The message was: ' e.message];
            logger.log(msg)
        end
    end
end

% 
% % save all open figures
% fig_handles = findall(0, 'Type', 'figure');
% if not(isempty(fig_handles))
%    for i = 1:length(fig_handles)
%        figname = num2str(fig_handles(i).Number);
%        try
%        savefig(fig_handles(i), fullfile(folder_path, 'plots', figname))
%        catch e
%         disp(['A problem occured while saving a fugure. Other data is saved. The message was: ' e.message])
%         logger.log(['a problem occured in Alkali_FID_measurement. The message was: ' e.message])
%        end
%    end
% end
