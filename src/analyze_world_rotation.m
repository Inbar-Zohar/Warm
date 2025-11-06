% ========================================================================
%   This script analyzes a world rotation measurement. It takes the time
%   trace of the rotation signal, calculates the Allan deviation and PSD,
%   and find all noise models that are consistent with both curves, plots
%   the results and saves them.
% ========================================================================

%% Defaults
% In principle the script should be run when a bunch of relevant fields are
% defined, but in case it's not, some defaults / backward compatability
if ~isfield(prm,'estimated_bw')
    prm.estimated_bw = 8;
end
if ~isfield(prm,'force_noise_term')
    prm.force_noise_term = [];
end
if ~isfield(prm,'default_tau_0')
    prm.default_tau_0 = 10;
end
if ~isfield(prm,'force_tau_0')
    prm.force_tau_0 = nan;
end
if ~isfield(prm,'max_t') || isnan(prm.max_t) || isempty(prm.max_t)
    prm.max_t = Inf;
end
if ~isfield(prm,'min_t') || isnan(prm.min_t) || isempty(prm.min_t)
    prm.min_t = -Inf;
end
if ~isfield(prm,'base_fig_num')
    prm.base_fig_num = 300;
end
if ~isfield(flg,'plot_individual_xenons')
    flg.plot_individual_xenons = 1;
end
if ~isfield(meas,'Omega_deg_hr') || (flg.plot_individual_xenons && ~isfield(meas,'Omega_129_deg_hr'))
    if isfield(meas,'Det131_PIDb_V') && isfield(meas,'Bz_V')
        %%
        [Omega_deg_hr,Bz,Omega_129_deg_hr,Omega_131_deg_hr] = ...
            get_Omega_and_Bz_from_signals(meas,prm);
        meas.Omega_deg_hr = Omega_deg_hr;
        meas.Bz = Bz;
        if flg.plot_individual_xenons
            meas.Omega_129_deg_hr = Omega_129_deg_hr;
            meas.Omega_131_deg_hr = Omega_131_deg_hr;
        end
    else
        error('No world rotation data');
    end
end
if ~isfield(meas,'fs')
    meas.fs = 1/diff(meas.time_s(1:2));
end

ga = (-7.44e7  )/(2*pi)/(1e4);              % [Hz/G] Xe129 gyromagnetic ratio
gb = (+2.2056e7)/(2*pi)/(1e4);              % [Hz/G] Xe131 gyromagnetic ratio
R = abs(ga / gb);


%% Analysis
% Extracting the Allan, PSD, and fitted noise terms
[tau_s, AD, AD_err, f, srPSD, srPSD_err, f_full, srPSD_full,...
    tau_range, f_range,...
    noise_term_vals, noise_term_errs, ...
    allan_noise_term_vals, allan_noise_term_errs, ...
    srPSD_noise_term_vals, srPSD_noise_term_errs] = ...
    get_noise_spectrum_and_parameters(...
    meas.Omega_deg_hr(meas.time_s > prm.min_t & meas.time_s < prm.max_t), meas.fs, prm.estimated_bw, prm.force_noise_term, 1);

% We plot and conveniently save one fitted noise model. By default it's the
% best fit, but we can force it to be any pre-selected set of parameters.
if ~isempty(prm.force_noise_term)
    ind2disp = find(ismember(isfinite(noise_term_vals(:,1:5)),prm.force_noise_term,'rows'));
else
    ind2disp = 1;
end
disp_noise_term_vals = noise_term_vals(ind2disp,:);
disp_noise_term_errs = noise_term_errs(ind2disp,:);
disp_res_string = iif(ind2disp > 1,'forced params','best fit');

% Get derived values: effective/noise bandwidth, effective bias
% instability, and optimal averaging time.
[BWeff,BWeff_err] = get_effective_BW(disp_noise_term_vals,disp_noise_term_errs);
[BIeff,BIeff_err,is_bound_BIeff,tau_min] = ...
    get_effective_BI(disp_noise_term_vals,disp_noise_term_errs,tau_range);
tau0 = iif(~is_bound_BIeff,tau_min,nan);


% Analyze individual Xenon signals
if flg.plot_individual_xenons
    %%
    ind_force_noise_terms = logical([1 1 1 1 1]);
    [~, AD_129, AD_err_129, ~, srPSD_129, srPSD_err_129, ~, ~,...
        ~, ~,...
        noise_term_vals_129, noise_term_errs_129] = ...
        get_noise_spectrum_and_parameters(...
        meas.Omega_129_deg_hr(meas.time_s > prm.min_t & meas.time_s < prm.max_t)/(1+R), meas.fs, ...
        prm.estimated_bw, ind_force_noise_terms, 1);
    if ~isempty(ind_force_noise_terms)
        ind2disp = find(ismember(isfinite(noise_term_vals_129(:,1:5)),ind_force_noise_terms,'rows'));
    else
        ind2disp = 1;
    end
    noise_term_vals_129 = noise_term_vals_129(ind2disp,:);
    noise_term_errs_129 = noise_term_errs_129(ind2disp,:);

    [~, AD_131, AD_err_131, ~, srPSD_131, srPSD_err_131, ~, ~,...
        ~, ~,...
        noise_term_vals_131, noise_term_errs_131] = ...
        get_noise_spectrum_and_parameters(...
        meas.Omega_131_deg_hr(meas.time_s > prm.min_t & meas.time_s < prm.max_t)*R/(1+R), meas.fs, ...
        prm.estimated_bw, ind_force_noise_terms, 1);
    if ~isempty(ind_force_noise_terms)
        ind2disp = find(ismember(isfinite(noise_term_vals_131(:,1:5)),ind_force_noise_terms,'rows'));
    else
        ind2disp = 1;
    end

    noise_term_vals_131 = noise_term_vals_131(ind2disp,:);
    noise_term_errs_131 = noise_term_errs_131(ind2disp,:);

    [this_PSD_full,f_full_Bz] = pwelch(meas.Bz - mean(meas.Bz), [], [], [], meas.fs);
    f_full_Bz = f_full_Bz(2:end-1);
    this_srPSD_full = sqrt(this_PSD_full(2:end-1));
    [f_Bz,srPSD_Bz] = logarithmic_smooth(f_full_Bz, this_srPSD_full, 'ieee-log');
end


%% Save stuff to meas structure
meas.tau_s_AD = tau_s;
meas.AD = AD;
meas.AD_err = AD_err;
meas.f_srPSD = f;
meas.srPSD = srPSD;
meas.srPSD_err = srPSD_err;
meas.tau_range = tau_range;
meas.f_range = f_range;

meas.AWN = disp_noise_term_vals(1);
meas.AWN_err = disp_noise_term_errs(1);
meas.ARW = disp_noise_term_vals(2);
meas.ARW_err = disp_noise_term_errs(2);
meas.BI = disp_noise_term_vals(3);
meas.BI_err = disp_noise_term_errs(3);
meas.RRW = disp_noise_term_vals(4);
meas.RRW_err = disp_noise_term_errs(4);
meas.RR = disp_noise_term_vals(5);
meas.RR_err = disp_noise_term_errs(5);
meas.AWNd = disp_noise_term_vals(6);
meas.AWNd_err = disp_noise_term_errs(6);

meas.BWeff = BWeff;
meas.BWeff_err = BWeff_err;
meas.BIeff = BIeff;
meas.BIeff_err = BIeff_err;
meas.tau_0 = tau0;

meas.noise_term_vals = noise_term_vals;
meas.noise_term_errs = noise_term_errs;
meas.allan_noise_term_vals = allan_noise_term_vals;
meas.allan_noise_term_errs = allan_noise_term_errs;
meas.srPSD_noise_term_vals = srPSD_noise_term_vals;
meas.srPSD_noise_term_errs = srPSD_noise_term_errs;


%% Plots
base_fig_num = prm.base_fig_num;

if 0
    %% Plot the comparison of the full and smoothed sqrt(PSD) curves
    this_fig = base_fig_num+10;
    meas = create_fig(meas,this_fig);
    plot_allan_srPSD_data(f_full, srPSD_full, 'srPSD', ...
        'plot_color', brighten(linecolor(1),0.5), 'lw', 0.5, 'legend_entry', 'Full');
    plot_allan_srPSD_data(f, srPSD, [], [], 'srPSD', ...
        'plot_color', linecolor(1), 'lw', 2, 'legend_entry', 'Downsampled');
    title({char(meas.DateTime_var), '\surdPSD smoothing'})
    set(gca,'XScale','log','Yscale','log');
    set(gcf,'Position', [    1180   200   560   350]);
end

if 1
    %% Plot one set of noise values (either the best fit or the preselected set)
    this_fig = base_fig_num+1;
    meas = create_fig(meas,this_fig);
    model_cmap = get_all_terms_cmap;
    % model_cmap = model_cmap(1,:);
    ax_single = plot_allan_srPSD_noise_wrapper(meas.tau_s_AD, meas.AD, meas.AD_err, ...
        meas.f_srPSD, meas.srPSD, meas.srPSD_err, ...
        meas.tau_range, meas.f_range,...
        [meas.AWN meas.ARW meas.BI meas.RRW meas.RR meas.AWNd], ...
        [meas.AWN_err meas.ARW_err meas.BI_err meas.RRW_err meas.RR meas.AWNd_err], ...
        'fig',this_fig, 'reset_fig',1, 'data_leg', [char(meas.DateTime_var) ', ' disp_res_string], 'plot_all', 1, ...
        'linecolor_data',linecolor(1), 'linecolor_model',model_cmap);
    BI_eff_range_for_plot = [max(BIeff/10,BIeff-BIeff_err) BIeff+BIeff_err] * sqrt(2*log(2)/pi);
    % plot(ax_single(1),tau0*[1 1],BI_eff_range_for_plot , '--k','HandleVisibility','off')
    BI_eff_for_plot = BIeff * sqrt(2*log(2)/pi);
    ax_single = get_ax_handles(this_fig);
    ax_single(1).XLim(1) = 0.1;
    ax_single(2).XLim(2) = 20;

    errorbar(ax_single(1),tau0,BIeff * sqrt(2*log(2)/pi), BIeff_err * sqrt(2*log(2)/pi),...
        'k','linewidth',2,'CapSize',3,'HandleVisibility','off');
    set(gcf,'Position', [    20   200   560   700]);
end

if 0
    %% Plot all noise sets of noise values that are consistent with both Allan and PSD
    this_fig = base_fig_num+2;
    meas = create_fig(meas,this_fig);
    ax_all = plot_allan_srPSD_noise_wrapper(tau_s, AD, AD_err, f, srPSD, srPSD_err, ...
        tau_range, f_range,...
        noise_term_vals, noise_term_errs,...
        'fig',this_fig, 'reset_fig',1, 'data_leg', 'Data', 'plot_all', 0, ...
        'show_title',0, 'show_legend',1);
    title(ax_all(1),[char(meas.DateTime_var) ', all fits'])
    set(gcf,'Position', [    600   200   560   700]);
end

tau0_for_plot = iif(isfinite(prm.force_tau_0),prm.force_tau_0, ...
    iif(isfinite(tau0),min(tau0, meas.time_s(end)/10),meas.time_s(end)/10));

corr_factor_less = 3;
smooth_param_corr_less = round(tau0_for_plot * meas.fs/corr_factor_less);
n_digits_less = max(2,ceil(log10(tau0_for_plot/corr_factor_less)));
label_smooth_less = sprintf('%s = %s s', '\tau_{avg}', num2str(tau0_for_plot/corr_factor_less,n_digits_less));
time_vec_smooth_less = meas.time_s(1:smooth_param_corr_less:end);
time_vec_smooth_less = time_vec_smooth_less(1:length(decimate(meas.time_s,smooth_param_corr_less)));
time_vec_smooth_less = time_vec_smooth_less + meas.time_s(floor(smooth_param_corr_less/2));

smooth_param_corr = round(tau0_for_plot * meas.fs);
n_digits = max(2,ceil(log10(tau0_for_plot)));
label_smooth = sprintf('%s = %s s', '\tau_{avg}', num2str(tau0_for_plot,n_digits));
time_vec_smooth = meas.time_s(1:smooth_param_corr:end);
time_vec_smooth = time_vec_smooth(1:length(decimate(meas.time_s,smooth_param_corr)));
time_vec_smooth = time_vec_smooth + meas.time_s(floor(smooth_param_corr/2));

analysis_range = [max(meas.time_s(1),prm.min_t) min(meas.time_s(end),prm.max_t)];


if 1
    %% Plot raw and smoothed Bz signal.

    this_fig = base_fig_num+4;
    meas = create_fig(meas,this_fig);
    try
        tiledlayout(2,1,'Padding','compact','TileSpacing','compact');
        ax(1) = nexttile;
        grid on; hold on; box on;
        ax(2) = nexttile;
        grid on; hold on; box on;
    catch
        for ind=1:2
            ax(ind) = subplot(2,1,ind);
            hold on; grid on; box on;
        end
    end
    axes(ax(1));
    plot(meas.time_s,meas.Bz*1e6,'DisplayName','Raw')
    if any(isfinite([prm.max_t prm.min_t]))
        yrange = ylim;
        ptch = patch('XData',analysis_range([1 1 2 2 1]),'YData',yrange([1 2 2 1 1]), ...
            'FaceColor','k','FaceAlpha',0.05,'EdgeColor','none');
        ax(1).Children = ax(1).Children([2:end 1]);
        ptch.HandleVisibility = 'off';
        ylim(yrange)
    end
    xlim(meas.time_s([1 end]))
    ylabel('dB_z (\muG)')
    title_str = 'dB_z smoothing';
    if ~is_bound_BIeff
        title_str = [title_str, ', \tau_{0} = ' num2str(tau0,max(2,ceil(log10(tau0)))) ' s'];
    end
    title({char(meas.DateTime_var), title_str})

    axes(ax(2));
    plot(nan,nan,'DisplayName','Raw')

    [time_vec_smooth_less, smoothed_Bz_less] = mymovmean(meas.time_s,meas.Bz,smooth_param_corr_less);
    plot(time_vec_smooth_less,smoothed_Bz_less*1e6,'LineWidth',1,'DisplayName',label_smooth_less);

    [time_vec_smooth, smoothed_Bz] = mymovmean(meas.time_s,meas.Bz,smooth_param_corr);
    plot(time_vec_smooth,smoothed_Bz*1e6,'LineWidth',2,'DisplayName',label_smooth);

    if any(isfinite([prm.max_t prm.min_t]))
        yrange = ylim;
        ptch = patch('XData',analysis_range([1 1 2 2 1]),'YData',yrange([1 2 2 1 1]), ...
            'FaceColor','k','FaceAlpha',0.05,'EdgeColor','none');
        ax(2).Children = ax(2).Children([2:end 1]);
        ptch.HandleVisibility = 'off';
        ylim(yrange)
    end
    xlim(meas.time_s([1 end]))
    xlabel('time (s)')
    ylabel('dB_z (\muG)')
    set(gcf,'Position', [    1180   200   560   700]);
    legend('location','best')
end
if 1
    %% Plot raw and smoothed Bz signal.

    this_fig = base_fig_num+8;
    meas = create_fig(meas,this_fig);
    try
        tiledlayout(2,1,'Padding','compact','TileSpacing','compact');
        ax(1) = nexttile;
        grid on; hold on; box on;
        ax(2) = nexttile;
        grid on; hold on; box on;
    catch
        for ind=1:2
            ax(ind) = subplot(2,1,ind);
            hold on; grid on; box on;
        end
    end
    axes(ax(1));
    plot(meas.time_s,meas.dESR*1e6,'DisplayName','Raw')
    if any(isfinite([prm.max_t prm.min_t]))
        yrange = ylim;
        ptch = patch('XData',analysis_range([1 1 2 2 1]),'YData',yrange([1 2 2 1 1]), ...
            'FaceColor','k','FaceAlpha',0.05,'EdgeColor','none');
        ax(1).Children = ax(1).Children([2:end 1]);
        ptch.HandleVisibility = 'off';
        ylim(yrange)
    end
    xlim(meas.time_s([1 end]))
    ylabel('dESR (\muG)')
    title_str = 'dESR smoothing';
    if ~is_bound_BIeff
        title_str = [title_str, ', \tau_{0} = ' num2str(tau0,max(2,ceil(log10(tau0)))) ' s'];
    end
    title({char(meas.DateTime_var), title_str})

    axes(ax(2));
    plot(nan,nan,'DisplayName','Raw')

    [time_vec_smooth_less, smoothed_ESR_less] = mymovmean(meas.time_s,meas.dESR,smooth_param_corr_less);
    plot(time_vec_smooth_less,smoothed_ESR_less*1e6,'LineWidth',1,'DisplayName',label_smooth_less);

    [time_vec_smooth, smoothed_ESR] = mymovmean(meas.time_s,meas.dESR,smooth_param_corr);
    plot(time_vec_smooth,smoothed_ESR*1e6,'LineWidth',2,'DisplayName',label_smooth);

    if any(isfinite([prm.max_t prm.min_t]))
        yrange = ylim;
        ptch = patch('XData',analysis_range([1 1 2 2 1]),'YData',yrange([1 2 2 1 1]), ...
            'FaceColor','k','FaceAlpha',0.05,'EdgeColor','none');
        ax(2).Children = ax(2).Children([2:end 1]);
        ptch.HandleVisibility = 'off';
        ylim(yrange)
    end
    xlim(meas.time_s([1 end]))
    xlabel('time (s)')
    ylabel('dESR (\muG)')
    set(gcf,'Position', [    1180   200   560   700]);
    legend('location','best')
end
if 1
    %% Plot raw and smoothed rotation signal.
    this_fig = base_fig_num+3;
    meas = create_fig(meas,this_fig);
    try
        tiledlayout(2,1,'Padding','compact','TileSpacing','compact');
        ax(1) = nexttile;
        grid on; hold on; box on;
        ax(2) = nexttile;
        grid on; hold on; box on;
    catch
        for ind=1:2
            ax(ind) = subplot(2,1,ind);
            hold on; grid on; box on;
        end
    end

    axes(ax(1));
    plot(meas.time_s,meas.Omega_deg_hr,'DisplayName','Raw')
    if any(isfinite([prm.max_t prm.min_t]))
        yrange = ylim;
        ptch = patch('XData',analysis_range([1 1 2 2 1]),'YData',yrange([1 2 2 1 1]), ...
            'FaceColor','k','FaceAlpha',0.05,'EdgeColor','none');
        ax(1).Children = ax(1).Children([2:end 1]);
        ptch.HandleVisibility = 'off';
        ylim(yrange)
    end
    xlim(meas.time_s([1 end]))
    ylabel('\Omega (deg/hour)')
    title_str = '\Omega smoothing';
    if ~is_bound_BIeff
        title_str = [title_str, ', \tau_{0} = ' num2str(tau0,max(2,ceil(log10(tau0)))) ' s'];
    end
    title({char(meas.DateTime_var), title_str})

    axes(ax(2));
    plot(nan,nan,'DisplayName','Raw')

    [time_vec_smooth_less, smoothed_Omega_less] = mymovmean(meas.time_s,meas.Omega_deg_hr,smooth_param_corr_less);
    plot(time_vec_smooth_less,smoothed_Omega_less,'LineWidth',1,'DisplayName',label_smooth_less)

    [time_vec_smooth, smoothed_Omega] = mymovmean(meas.time_s,meas.Omega_deg_hr,smooth_param_corr);
    plot(time_vec_smooth,smoothed_Omega,'LineWidth',2,'DisplayName',label_smooth)

    if any(isfinite([prm.max_t prm.min_t]))
        yrange = ylim;
        ptch = patch('XData',analysis_range([1 1 2 2 1]),'YData',yrange([1 2 2 1 1]), ...
            'FaceColor','k','FaceAlpha',0.05,'EdgeColor','none');
        ax(2).Children = ax(2).Children([2:end 1]);
        ptch.HandleVisibility = 'off';
        ylim(yrange)
    end

    xlim(meas.time_s([1 end]))
    xlabel('time (s)')
    ylabel('\Omega (deg/hour)')
    set(gcf,'Position', [    1180   200   560   700]);
    legend('location','best')
end


if flg.plot_individual_xenons && isfield(meas,'Omega_129_deg_hr')
    %%
    if 1
        %%
        this_fig = base_fig_num+6;
        meas = create_fig(meas,this_fig);
        t = tiledlayout(1,1,"TileSpacing","compact","Padding","compact");
        nexttile;
        set(gca,'XScale','log','YScale','log'); hold on; grid on;
        plot_allan_srPSD_data(f,srPSD_129,'srPSD','lw',2,'legend_entry','\Omega_{129} / (1+R)','plot_color',linecolor(1));
        plot_allan_srPSD_data(f,srPSD_131,'srPSD','lw',2,'legend_entry','\Omega_{131} \cdot R/(1+R)','plot_color',linecolor(2));
        plot_allan_srPSD_data(f,srPSD,'srPSD','lw',2,'legend_entry','\Omega_{rot}','plot_color',linecolor(3));
        plot_allan_srPSD_data(f_Bz,srPSD_Bz * abs(ga)/(1+R)*360*3600/60,'srPSD','lw',2,'legend_entry','\Omega_{B}','plot_color',color_vec('k'));

        if isfield(prm,'WP'); plot_theory_psd_tmp; end

        xl = extend_range(f(1),f(end),1);
        xl(2) = 30;
        xlim(xl)

        yl1 = extend_range(srPSD_129_calc(1)/(1+R)/5,max(srPSD_129)*60,1);
        yl2 = extend_range(srPSD_131_calc(1)*R/(1+R)/5,max(srPSD_131)*60,1);
        yl = [min([yl1 yl2]) max([yl1 yl2])];
        ylim(yl)

        % legend('location','best')

        add_linked_secondary_yaxis(gca,@(x) x / (abs(ga)/(1+R)*360*3600));
        yyaxis right;
        ylabel('\surdPSD_{B_z} (G/\surdHz)')
        set(gca,'YScale','log')
        yyaxis left;

        title({char(meas.DateTime_var),'Normalized Xenons and \Omega/B_z',theory_title_str})

        set(gcf,'Position', [    600   200   560   600]);
    end
    if 1
        %%
        this_fig = base_fig_num+7;
        meas = create_fig(meas,this_fig);
        t = tiledlayout(2,2,"TileSpacing","compact","Padding","compact");
        nexttile; set(gca,'XScale','log','YScale','log'); grid on;
        all_ax(1,1) = gca;
        nexttile; set(gca,'XScale','log','YScale','log'); grid on;
        all_ax(1,2) = gca;
        nexttile; set(gca,'XScale','log','YScale','log'); grid on;
        all_ax(2,1) = gca;
        nexttile; set(gca,'XScale','log','YScale','log'); grid on;
        all_ax(2,2) = gca;
        t.Title.String = char(meas.DateTime_var);
        t.Title.FontWeight = 'bold';
        ax_129 = plot_allan_srPSD_noise_wrapper(tau_s, AD_129, AD_err_129, f, srPSD_129, srPSD_err_129, ...
            tau_range, f_range,...
            noise_term_vals_129, noise_term_errs_129,...
            'fig',all_ax(:,1), 'reset_fig',0, 'data_leg', 'Xe129 / (1+R)', 'plot_all', 1, ...
            'linecolor_data',linecolor(1), 'linecolor_model',model_cmap);

        % xlim(ax_129(1),[3e-2 3e2])
        % ylim(ax_129(1),[3e0 3e2])
        % xlim(ax_129(2),[1e-2 3e1])
        % ylim(ax_129(2),[3e0 10e2])

        ax_131 = plot_allan_srPSD_noise_wrapper(tau_s, AD_131, AD_err_131, f, srPSD_131, srPSD_err_131, ...
            tau_range, f_range,...
            noise_term_vals_131, noise_term_errs_131,...
            'fig',all_ax(:,2), 'reset_fig',0, 'data_leg', 'Xe131 \cdot R/(1+R)', 'plot_all', 1, ...
            'linecolor_data',linecolor(1), 'linecolor_model',model_cmap);


        % xlim(ax_131(1),ax_129(1).XLim)
        % ylim(ax_131(1),ax_129(1).YLim)
        % xlim(ax_131(2),ax_129(2).XLim)
        % ylim(ax_131(2),ax_129(2).YLim)
        set(gcf,'Position', [    20   200   560*20/10   700]);
    end
    if 0
        %%
        this_fig = base_fig_num+8;
        meas = create_fig(meas,this_fig);
        set(gca,'XScale','log','YScale','log')
        plot_allan_srPSD_data(f,srPSD_129*(1+R),'srPSD','lw',2,'legend_entry','\Omega_{129}','plot_color',linecolor(1));
        plot_allan_srPSD_data(f,srPSD_131*(1+R)/R,'srPSD','lw',2,'legend_entry','\Omega_{131}','plot_color',linecolor(2));
        xlim([1e-2 1e2])
        % ylim([3e1 3e3])
        title({char(meas.DateTime_var),'Raw Xenons'})
    end
end
