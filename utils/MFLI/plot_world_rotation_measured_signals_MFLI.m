
[srPSD_Xe_129_Theta,f] = pwelch_and_smooth(meas.Xe_129_Theta,meas.time_s);
srPSD_Xe_131_Theta = pwelch_and_smooth(meas.Xe_131_Theta,meas.time_s);

srPSD_Xe_129_R = pwelch_and_smooth(meas.Xe_129_R,meas.time_s);
srPSD_Xe_131_R = pwelch_and_smooth(meas.Xe_131_R,meas.time_s);

srPSD_Xe_129_Y = pwelch_and_smooth(meas.Xe_129_Y,meas.time_s);
srPSD_Xe_131_Y = pwelch_and_smooth(meas.Xe_131_Y,meas.time_s);

srPSD_Xe_129_Omega_deg_hr = pwelch_and_smooth(meas.Omega_129_deg_hr,meas.time_s);
srPSD_Xe_131_Omega_deg_hr = pwelch_and_smooth(meas.Omega_131_deg_hr,meas.time_s);

plot_sim = 0;
if isfield(meas,'Sim_Theta')
    plot_sim = 1;
    srPSD_Sim_Theta = pwelch_and_smooth(meas.Sim_Theta,meas.time_s);
    srPSD_Sim_R = pwelch_and_smooth(meas.Sim_R,meas.time_s);
    srPSD_Sim_Y = pwelch_and_smooth(meas.Sim_Y,meas.time_s);
end
%%
figure(105); clf;
t=tiledlayout(4,2,"TileSpacing","compact","Padding","compact");
title(t,{char(meas.DateTime_var), 'All time traces'})

time_unit_min = 0;
factor = iif(time_unit_min,60,1);
xlbl = sprintf('time (%s)',iif(time_unit_min,'min','s'));
xlabel(t,xlbl)
movmean_param = meas.fs*3;


nexttile; hold all; grid on;
plot(meas.time_s/factor, meas.Xe_129_Y - mean(meas.Xe_129_Y), 'color',linecolor(1))
plot(meas.time_s/factor, meas.Xe_131_Y - mean(meas.Xe_131_Y), 'color',linecolor(2))
if plot_sim; plot(meas.time_s/factor, meas.Sim_Y - mean(meas.Sim_Y), 'color',linecolor(3)); end
xlim(meas.time_s([1 end])/factor)
ylabel('Y (V)')

nexttile; hold all; grid on;
plot(meas.time_s/factor, movmean(meas.Xe_129_Y,movmean_param) - mean(meas.Xe_129_Y), 'color',brighten(linecolor(1),-0.3));
plot(meas.time_s/factor, movmean(meas.Xe_131_Y,movmean_param) - mean(meas.Xe_131_Y), 'color',brighten(linecolor(2),-0.3));
plot(meas.time_s/factor, movmean(meas.Xe_131_Y,movmean_param) - mean(meas.Xe_131_Y), 'color',brighten(linecolor(2),-0.3));
if plot_sim; plot(meas.time_s/factor, movmean(meas.Sim_Y,movmean_param) - mean(meas.Sim_Y), 'color',brighten(linecolor(3),-0.3));end

xlim(meas.time_s([1 end])/factor)

nexttile; hold all; grid on;
plot(meas.time_s/factor, meas.Xe_129_R ./ mean(meas.Xe_129_R) - 1, 'color',linecolor(1))
plot(meas.time_s/factor, meas.Xe_131_R ./ mean(meas.Xe_131_R) - 1, 'color',linecolor(2))
if plot_sim; plot(meas.time_s/factor, meas.Sim_R ./ mean(meas.Sim_R) - 1, 'color',linecolor(3)); end
xlim(meas.time_s([1 end])/factor)
ylabel('dR/R')

nexttile; hold all; grid on;
plot(meas.time_s/factor, movmean(meas.Xe_129_R,movmean_param) ./ mean(meas.Xe_129_R) - 1, 'color',brighten(linecolor(1),-0.3));
plot(meas.time_s/factor, movmean(meas.Xe_131_R,movmean_param) ./ mean(meas.Xe_131_R) - 1, 'color',brighten(linecolor(2),-0.3));
if plot_sim; plot(meas.time_s/factor, movmean(meas.Sim_R,movmean_param) ./ mean(meas.Sim_R) - 1, 'color',brighten(linecolor(3),-0.3));end
xlim(meas.time_s([1 end])/factor)

nexttile; hold all; grid on;
plot(meas.time_s/factor, meas.Xe_129_Theta - mean(meas.Xe_129_Theta), 'color',linecolor(1))
plot(meas.time_s/factor, meas.Xe_131_Theta - mean(meas.Xe_131_Theta), 'color',linecolor(2))
if plot_sim; plot(meas.time_s/factor, meas.Sim_Theta - mean(meas.Sim_Theta), 'color',linecolor(3)); end
xlim(meas.time_s([1 end])/factor)
ylabel('\Theta (deg)')

nexttile; hold all; grid on;
plot(meas.time_s/factor, movmean(meas.Xe_129_Theta,movmean_param) - mean(meas.Xe_129_Theta), 'color',brighten(linecolor(1),-0.3));
plot(meas.time_s/factor, movmean(meas.Xe_131_Theta,movmean_param) - mean(meas.Xe_131_Theta), 'color',brighten(linecolor(2),-0.3));
if plot_sim; plot(meas.time_s/factor, movmean(meas.Sim_Theta,movmean_param) - mean(meas.Sim_Theta), 'color',brighten(linecolor(3),-0.3));end
xlim(meas.time_s([1 end])/factor)

nexttile; hold all; grid on;
plot(meas.time_s/factor, meas.Omega_129_deg_hr - mean(meas.Omega_129_deg_hr), 'color',linecolor(1))
plot(meas.time_s/factor, meas.Omega_131_deg_hr - mean(meas.Omega_131_deg_hr), 'color',linecolor(2))
xlim(meas.time_s([1 end])/factor)
ylabel('\Omega (deg/hr)')
xlabel(xlbl)

nexttile; hold all; grid on;
plot(meas.time_s/factor, movmean(meas.Omega_129_deg_hr,movmean_param) - mean(meas.Omega_129_deg_hr), 'color',brighten(linecolor(1),-0.3));
plot(meas.time_s/factor, movmean(meas.Omega_131_deg_hr,movmean_param) - mean(meas.Omega_131_deg_hr), 'color',brighten(linecolor(2),-0.3));
xlim(meas.time_s([1 end])/factor)


%%
figure(106); clf; %set(gca,'XScale','log','Yscale','log'); hold all; grid on;
% this_fig = base_fig_num+1002;
% meas = create_fig(meas,this_fig);

fit_range = [1 5];
isotopes = {'Xe_129','Xe_131'}; if plot_sim; isotopes{end+1} = 'Sim'; end
params = {'Y','R','Theta'};
fit_vals = nan(length(isotopes),length(params));
for ind_i = 1:length(isotopes)
    iso = isotopes{ind_i};
    for ind_p = 1:length(params)
        par = params{ind_p};
        eval(sprintf('fit_vals(ind_i,ind_p) = mean(srPSD_%s_%s(f > fit_range(1) & f<fit_range(2)));',...
            iso,par));
    end
end
fit_vals = fit_vals .* repmat([1e6 1e6 1e3],length(isotopes),1) / sqrt(2);
round(fit_vals,1)
fig = gcf;
fig.Position([4]) = [740];
if fig.Position(2) > 400; fig.Position(2) = 200; end
 
t=tiledlayout(3,1,"TileSpacing","compact","Padding","compact");
title(t,{char(meas.DateTime_var), 'All PSD'})
 
xl = [5e-3 6e1];
 
nexttile; set(gca,'XScale','log','Yscale','log'); hold all; grid on;

plot(f,srPSD_Xe_129_Y,'-','LineWidth',2,'DisplayName','Y_{129}','Color',linecolor(1))
plot(f,srPSD_Xe_131_Y,'-' ,'LineWidth',2,'DisplayName','Y_{131}','Color',linecolor(2))
if plot_sim; plot(f,srPSD_Sim_Y,'-' ,'LineWidth',2,'DisplayName','Y_{Sim}','Color',linecolor(3)); end

plot(f,srPSD_Xe_129_R,'--','LineWidth',2,'DisplayName','R_{129}','Color',brighten(linecolor(1),-0.3))
plot(f,srPSD_Xe_131_R,'--' ,'LineWidth',2,'DisplayName','R_{131}','Color',brighten(linecolor(2),-0.3))
if plot_sim; plot(f,srPSD_Sim_R,'--' ,'LineWidth',2,'DisplayName','R_{Sim}','Color',brighten(linecolor(3),-0.3)); end

% 
ylabel('\surdPSD_{sig} (V/\surdHz)')
legend('Location','best')
xlim(xl)
% % ylim([1e-8 5e-5])
% 
nexttile; set(gca,'XScale','log','Yscale','log'); hold all; grid on;
plot(f,srPSD_Xe_129_Theta ,'LineWidth',2,'DisplayName','\Theta_{129}','Color',linecolor(1))
plot(f,srPSD_Xe_131_Theta ,'LineWidth',2,'DisplayName','\Theta_{131}','Color',linecolor(2))
if plot_sim; plot(f,srPSD_Sim_Theta ,'LineWidth',2,'DisplayName','\Theta_{Sim}','Color',linecolor(3)); end

ylabel('\surdPSD_{\Theta} (deg/\surdHz)')
legend('Location','best')
xlim(xl)
% ylim([1e-4 1e0])

nexttile; set(gca,'XScale','log','Yscale','log'); hold all; grid on;
plot(f,srPSD_Xe_129_Omega_deg_hr ,'LineWidth',2,'DisplayName','\Omega_{129}','Color',linecolor(1))
plot(f,srPSD_Xe_131_Omega_deg_hr ,'LineWidth',2,'DisplayName','\Omega_{131}','Color',linecolor(2))
% plot(f,srPSD_PID_ESR_V_srHz ,'LineWidth',2,'DisplayName','PID_{ESR}','Color',linecolor(3))
% 
ylabel('\surdPSD_{\Omega} (deg/hr/\surdHz)')
legend('Location','best')
xlim(xl)
% 

xlabel('f (Hz)')


        % yrange = ylim;
        % if any(isfinite([prm.max_t prm.min_t]))
        %     ptch = patch('XData',analysis_range([1 1 2 2 1]),'YData',yrange([1 2 2 1 1]), ...
        %         'FaceColor','k','FaceAlpha',0.05,'EdgeColor','none');
        %     this_ax = gca;
        %     this_ax.Children = this_ax.Children([2:end 1]);
        %     ptch.HandleVisibility = 'off';
        % end
        % ylim(yrange)
