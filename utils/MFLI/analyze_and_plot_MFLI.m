base_folder = 'C:\Users\Constantine\Documents\Zurich Instruments\LabOne\WebServer';

session = 'session_20250821_104720_23';
stream = 'stream_002';


%% 
save_path = 'D:\Synced team folder\Data Archive\DorB\All Data\Currently synced data';
main_path = 'D:\Synced team folder\SynologyDrive\Dor B\DorB_operational';

%% Loop data
substitute_if_does_not_exist('is_part_of_loop',0);

%% Experiment description
exp_description = ['']; % 

%%
gongData = load('gong.mat');
gongObj = audioplayer(gongData.y,gongData.Fs);



%% 
[meas, prm, flg] = extract_raw_data_from_MFLI_folder2(fullfile(base_folder,session,stream), []); 
plot_world_rotation_measured_signals_MFLI; 

%%
flg.save_all_open_figs = 0; % when 0 saves only figures that are created using create_fig
flg.is_part_of_loop = is_part_of_loop;

% save data
flg.save_data = 1; 

%% Init Logger
logger = ExperimentLogger(fullfile(main_path, 'log'),[],~flg.is_part_of_loop);
if 1
    logger.log('#########################################################')
    logger.log('Experiment Description:')
    logger.log(exp_description)
    if is_part_of_loop
        logger.log(loop_parse_loop_data(prm.loop_data));
    end

    logger.log('#########################################################')
end


%%
if 1
    %%
    % prm.min_t = 2000;
    prm.max_t = 3000;
    % prm.force_noise_term = logical([1 1 1 0 0]);
    % 
    % prm.min_t = 2.5e4;
    % prm.max_t = 7000;
    % prm.force_noise_term = logical([1 1 0 0 1]);

    % prm.min_t = 7200;
    % prm.force_noise_term = logical([1 1 0 1 0]);

    % prm.force_noise_term = logical([1 0 1 0 1]);
    prm.estimated_bw = 8;
    % prm.base_fig_num = 600;
    analyze_world_rotation
    if 0
        %%
        dR = meas.Sim_R ./ mean(meas.Sim_R) - 1;
        inds = meas.time_s > prm.min_t & meas.time_s > prm.max_t;
        myfig(prm.base_fig_num + 50);
        myscatter(movmean(dR(inds),movmean_param)*10/5e-3, ...
            movmean(meas.Omega_deg_hr(inds),smooth_param_corr*3),[],1:length(dR),'alpha',0.3,'ms',3);

        % meas.Omega_deg_hr = meas.Omega_deg_hr  - 860*dR
        % prm.base_fig_num = 400;
        % analyze_world_rotation
    end
end

if 0
    %%
    load(fullfile(base_folder,session,stream,'stream_00000.mat'))

end

%% Save Data
logger.log('Done')

if flg.save_data
    % save_current_working_point
    save_everything
    logger.log('Saved everything')
end

%% reset logger

reset_all = 0;  % set reset to zero
if ~flg.is_part_of_loop || prm.loop_data.ind_run == prm.loop_data.n_runs_total
    playblocking(gongObj);
end
logger.close_log