% ========================================================================
%   This script is for measuring Xenon FID, computing the T2s and resonance 
%   frequncies from the data. The code assumes Bz, Bx and By are set, 
%   NMRs are working and the PIDs are disabled. 
%   The signal is measured from the scope's channel 2 which is
%   assumed to be the signal after the ESR LIA.
%
% ========================================================================
%
%   TODO:
%   1. Get magnetization estimations from a fit rather than the 
%      "numeric lock-in"
%   2. 
% ========================================================================


% logger.log('Xenon FID measurement initialized.')
% if ~flg.is_part_of_loop
%     disp('[+] Starting a Xenon FID measurement');
% end

Vscale_Xe_FID = 0.15;

scope2.Reset
scope2.OperationComplete;
scope2.setChDisp(2);
scope2.setVscale(2,Vscale_Xe_FID);
scope2.setVoffset(2,0);

if flg.measure_magnetometer_sensitivity
    % apply sine modulation to By far from NMR resonance
    DorB.By.apply('Sin', 220, 5e-3, 0, prm.WP.fields.By0);
%     AG2.Sin(2, 220, 5e-3, 0, prm.WP.fields.By0);
    [~,wvfm] = DorB.By.apply('Output',[],true);
%     [~, wvfm] = AG2.Output(2,[],true);
    prm.WP.Vy_modulation.frequency = wvfm{1}.p(1);
    prm.WP.Vy_modulation.Amplitude = wvfm{1}.p(2);% / 2;
    prm.WP.Vy_modulation.Bamp = prm.WP.Vy_modulation.Amplitude * prm.WP.coils.dBy_dVy / 2;
else
    prm.WP.Vy_modulation.frequency = 0;
    prm.WP.Vy_modulation.Bamp = 0;
end

% set scope or DAQ for Xenon FID measurement
if flg.use_DAQ_4_Xe_FID
    error('using DAQ for Xe FID is not supported yet')
else
    Tscale_Xe_FID = 5;
    Tdelay = 4.0 * Tscale_Xe_FID;
    scope2.setTscale(Tscale_Xe_FID);
    scope2.setTdelay(Tdelay);
    scope2.setTref('CENTer');
    scope2.TrigSource('EXT');
    scope2.Single;
    pause(5 * Tscale_Xe_FID);
    
    % start measurement
    scope2.ForceTrig;
end

%% Start experiment 
% kill NMRs
AG1.OutputOFF(1); 
AG2.OutputOFF(1);
pause(5 * Tscale_Xe_FID);
scope2.readyToRead;

% start NMRs 
AG1.OutputON(1); 
AG2.OutputON(1);

if flg.measure_magnetometer_sensitivity
% kill the By sine modulation
AG2.DC(2, prm.WP.fields.By0)
end

%% read data
% read signal from scope
[t, v] = scope2.Read(2);
% t = t(1:round(length(t) / 10000):end);
% v = v(1:round(length(v) / 10000):end);

% release scope
scope2.setTref('LEFT'); 
scope2.setTdelay(0);
scope2.TrigSource('LINE'); 
scope2.Run;


%% Save data
meas.FID_Xe_t = t;
meas.FID_Xe_v = v(1,:);
clear t v

logger.log('Xenon FID data processing.')
% if ~flg.is_part_of_loop
%     disp('[+] Xenon FID data processing');
% end

%% Analyze
% set local variables


t = meas.FID_Xe_t;
Magnetometer = meas.FID_Xe_v;

% estimating amplitudes for coarse fitting
if flg.measure_magnetometer_sensitivity
prm.WP.sensitivity = abs(2 * trapz(t, Magnetometer .* exp(-2i * pi *...
    prm.WP.Vy_modulation.frequency * t))  / (t(end) - t(1))) /...
    prm.WP.Vy_modulation.Bamp;
else
    prm.WP.sensitivity = 1;
end
prm.WP.Xe129.B = abs(2 * trapz(t(t < 0), Magnetometer(t < 0) .* exp(-2i *...
    pi * prm.WP.Xe129.drive_nom.p(1) * t(t < 0))) / (max(t(t < 0)) -...
    min(t(t < 0))) / prm.WP.sensitivity);
prm.WP.Xe131.B = abs(2 * trapz(t(t < 0), Magnetometer(t < 0) .* exp(-2i *...
    pi * prm.WP.Xe131.drive_nom.p(1) * t(t < 0))) / (max(t(t < 0)) -...
    min(t(t < 0))) / prm.WP.sensitivity);
sine3model = fittype([...
    num2str(prm.WP.Xe129.B * prm.WP.sensitivity) '*exp(-G129*x)*cos(2*pi*f129*x+p129) + ' ...
    num2str(prm.WP.Xe131.B * prm.WP.sensitivity) '*exp(-G131*x)*cos(2*pi*f131*x+p131) + ' ...
    num2str(prm.WP.Vy_modulation.Bamp * prm.WP.sensitivity) '*cos(2*pi*' num2str(prm.WP.Vy_modulation.frequency) '*x+py) + dc'],...
    'coefficients',{'G129','f129','p129','G131','f131','p131','py','dc'});


% A time buffer for when the signal start to decay
t0_for_fit = 5e-3;

% coarse fit using amplitudes estimation
sine3fit = fit(t(t > t0_for_fit).', Magnetometer(t > t0_for_fit).', sine3model, 'StartPoint',...
    [0.12, prm.WP.Xe129.drive_nom.p(1), 0.5, 0.12, prm.WP.Xe131.drive_nom.p(1), 0.5, 0.5, mean(Magnetometer)]);
sine3model = fittype([...
    'A129*exp(-G129*x)*cos(2*pi*f129*x+p129) + ' ...
    'A131*exp(-G131*x)*cos(2*pi*f131*x+p131) + ' ...
    'Ay*cos(2*pi*' num2str(prm.WP.Vy_modulation.frequency) '*x+py) + dc'],...
    'coefficients',{'A129','G129','f129','p129','A131','G131','f131','p131','Ay','py','dc'});

% fine fit using coarse fit results
[sine3fit,gof] = fit(t(t > t0_for_fit).', Magnetometer(t > t0_for_fit).', sine3model, 'StartPoint',...
    [prm.WP.Xe129.B * prm.WP.sensitivity, sine3fit.G129, sine3fit.f129, sine3fit.p129,...
     prm.WP.Xe131.B * prm.WP.sensitivity, sine3fit.G131, sine3fit.f131, sine3fit.p131,...
     prm.WP.Vy_modulation.Bamp * prm.WP.sensitivity, sine3fit.py, sine3fit.dc]);

 fitresult_cint = confint(sine3fit, erf(1 / sqrt(2)));
 fitresult_cint = diff(fitresult_cint,1)/2;

 
% save relevant parameters 
prm.WP.Xe129.G2 = sine3fit.G129; 
prm.WP.Xe129.G2err = fitresult_cint(2); 
prm.WP.Xe129.fFID = sine3fit.f129;
prm.WP.Xe131.G2 = sine3fit.G131; 
prm.WP.Xe131.G2err = fitresult_cint(6); 
prm.WP.Xe131.fFID = sine3fit.f131;

T2_Xe129 = 1 / prm.WP.Xe129.G2;
T2_Xe129_err = (prm.WP.Xe129.G2err / prm.WP.Xe129.G2) * T2_Xe129;
T2_Xe131 = 1 / prm.WP.Xe131.G2;
T2_Xe131_err = (prm.WP.Xe131.G2err / prm.WP.Xe131.G2) * T2_Xe131;



% plot Xenon FID fit
cmap = lines;

typical_data_dt = 0.25e-3;
downsample_factor_data = ceil(typical_data_dt/diff(t(1:2)));
fit_to_data_factor = 2;

inds_disp_data = 1:downsample_factor_data:length(t);
t_disp_data = t(inds_disp_data);
Magnetometer_disp = Magnetometer(inds_disp_data);

t_disp_fit = linspace(t_disp_data(1),t_disp_data(end),length(t_disp_data)*5);
t_disp_fit = t_disp_fit(t_disp_fit > t0_for_fit);
fit_disp = feval(sine3fit, t_disp_fit);
fit_envelope_fixed = sine3fit.dc + [1;-1]*(sine3fit.Ay + 0*t_disp_fit) ;
fit_envelope_129 = sine3fit.A129 *exp(-sine3fit.G129*t_disp_fit);
fit_envelope_131 = sine3fit.A131 *exp(-sine3fit.G131*t_disp_fit);

fit_envelope_129_disp = fit_envelope_fixed + [1;-1]*fit_envelope_129;
fit_envelope_131_disp = fit_envelope_fixed + [1;-1]*fit_envelope_131;
fit_envelope_tot_disp = fit_envelope_fixed + [1;-1]*(fit_envelope_129+fit_envelope_131);

yl = sine3fit.dc + 1.1 * [-1 1]*(sine3fit.Ay + sine3fit.A129 +  sine3fit.A131);
yl(1) = min(yl(1), min(Magnetometer) - 0.1 * (max(Magnetometer) - min(Magnetometer)));
yl(2) = max(yl(2), max(Magnetometer) + 0.1 * (max(Magnetometer) - min(Magnetometer)));

zoom_t_range = 20e-3;
zoom_t_start = [t0_for_fit  round(t(end)/4) round(t(end)/2) round(0.9*t(end))];
zoom_inds_data = [];
zoom_inds_fit = [];

for ind=1:length(zoom_t_start)
    this_inds = find(t_disp_data >=  zoom_t_start(ind) & t_disp_data <=  zoom_t_start(ind) + zoom_t_range);
    zoom_inds_data(ind,:) = this_inds([1 end]);
    this_inds = find(t_disp_fit >=  zoom_t_start(ind) & t_disp_fit <=  zoom_t_start(ind) + zoom_t_range);
    zoom_inds_fit(ind,:) = this_inds([1 end]);
end

meas = create_fig(meas,5);
set(gcf,'position', [300 300 900 650]);

subplot(2,length(zoom_t_start),1:length(zoom_t_start))
hold all
grid on;
pos = get(gca,'position');
set(gca,'position', pos + [0 -0.09 0 0.02]);

for ind=1:length(zoom_t_start)
    patch('Xdata',zoom_t_start(ind)+zoom_t_range*[0 0 1 1], 'Ydata',yl([1 2 2 1]), ...
        'facecolor', 0.5*ones(1,3), 'facealpha', 0.3, 'edgecolor', 0.5*ones(1,3)) 

end

plot(t_disp_data, Magnetometer_disp,'.','color',cmap(1,:))
plot(t_disp_fit, fit_disp, '-','color',cmap(2,:));
plot(t_disp_fit, fit_envelope_fixed, ':','color','k','linewidth', 2);
plot(t_disp_fit, fit_envelope_129_disp, ':','color',cmap(3,:),'linewidth', 2);
plot(t_disp_fit, fit_envelope_131_disp, ':','color',cmap(4,:),'linewidth', 2);
plot(t_disp_fit, fit_envelope_tot_disp, '--','color','k','linewidth', 2);


if flg.measure_magnetometer_sensitivity
    magsen_str = ['magnetometer sensitivity = ' num2str(prm.WP.sensitivity / 1e3) ' V/mG\n'];
else
    magsen_str = 'magnetometer sensitivity = NA\n';
end

ttl_str = sprintf([char(meas.DateTime_var), '\n', magsen_str,...
'B[^{129}Xe]=' num2str(prm.WP.Xe129.B * 1e6, 3) ' \\muG , T_2[^{129}Xe]=' numerr2str(T2_Xe129,T2_Xe129_err) ' s , '...
'B[^{131}Xe]=' num2str(prm.WP.Xe131.B * 1e6, 3) ' \\muG , T_2[^{131}Xe]=' numerr2str(T2_Xe131,T2_Xe131_err) ' s\n'...
'f[^{129}Xe]=' num2str(prm.WP.Xe129.drive_nom.p(1)) ' Hz , f[^{131}Xe]=' num2str(prm.WP.Xe131.drive_nom.p(1)) ' Hz']);
title(ttl_str); 
xlabel('time (s)'); 
ylabel('signal (V)'); 
xlim(t([1 end]));
ylim(yl);

yyaxis right;
xlim(t([1 end]));
set(gca,'YColor','black');
ylabel('signal (\muG)'); 
ylim((yl-sine3fit.dc) /prm.WP.sensitivity *1e6);


for ind=1:length(zoom_t_start)
    subplot(2,length(zoom_t_start),length(zoom_t_start)+ind)
    hold all
    grid on;
    pos = get(gca,'position');
    set(gca,'position', pos - [0 0.04 0 0]);

    plot((t_disp_data(zoom_inds_data(ind,1):zoom_inds_data(ind,2))-zoom_t_start(ind))*1e3, ...
        Magnetometer_disp(zoom_inds_data(ind,1):zoom_inds_data(ind,2)),...
        '.','color',cmap(1,:),'markersize',8)
    plot((t_disp_fit(zoom_inds_fit(ind,1):zoom_inds_fit(ind,2))-zoom_t_start(ind))*1e3, ...
        fit_disp(zoom_inds_fit(ind,1):zoom_inds_fit(ind,2)), ...
        '-','color',cmap(2,:));
    plot((t_disp_fit(zoom_inds_fit(ind,1):zoom_inds_fit(ind,2))-zoom_t_start(ind))*1e3, ...
        fit_envelope_fixed(:,zoom_inds_fit(ind,1):zoom_inds_fit(ind,2)), ...
        ':','color','k','linewidth', 2);    
    plot((t_disp_fit(zoom_inds_fit(ind,1):zoom_inds_fit(ind,2))-zoom_t_start(ind))*1e3, ...
        fit_envelope_129_disp(:,zoom_inds_fit(ind,1):zoom_inds_fit(ind,2)), ...
        ':','color',cmap(3,:),'linewidth', 2);
    plot((t_disp_fit(zoom_inds_fit(ind,1):zoom_inds_fit(ind,2))-zoom_t_start(ind))*1e3, ...
        fit_envelope_131_disp(:,zoom_inds_fit(ind,1):zoom_inds_fit(ind,2)), ...
        ':','color',cmap(4,:),'linewidth', 2);
    plot((t_disp_fit(zoom_inds_fit(ind,1):zoom_inds_fit(ind,2))-zoom_t_start(ind))*1e3, ...
        fit_envelope_tot_disp(:,zoom_inds_fit(ind,1):zoom_inds_fit(ind,2)), ...
        '--','color','k','linewidth', 2);

    xlim([0, zoom_t_range*1e3]);
    ylim(yl);

    if zoom_t_start(ind) >= 0.1
        xlabel(['t-' num2str(zoom_t_start(ind)) 's (ms)']); 
    else
        xlabel(['t-' num2str(zoom_t_start(ind)*1e3 ) 'ms (ms)']); 
    end
    if ind == 1
        ylabel('signal (V)');
    else
%         set(gca,'YTickLabel',[])
    end
    if ind == length(zoom_t_start)
        yyaxis right;
        xlim([0, zoom_t_range*1e3]);
        set(gca,'YColor','black');
        ylabel('signal (\muG)');
        ylim((yl-sine3fit.dc) /prm.WP.sensitivity *1e6);
        
    end
   
end
%%

% 'A129','G129','f129','p129','A131','G131','f131','p131','Ay','py','dc'
% prm.WP.Xe129.G2 = sine3fit.G129; 
% prm.WP.Xe129.G2err = fitresult_cint(2); 
% prm.WP.Xe129.fFID = sine3fit.f129;
% prm.WP.Xe131.G2 = sine3fit.G131; 
% prm.WP.Xe131.G2err = fitresult_cint(6); 
% prm.WP.Xe131.fFID = sine3fit.f131;

% T2_Xe129 = 1 / prm.WP.Xe129.G2;
% T2_Xe129_err = (prm.WP.Xe129.G2err / prm.WP.Xe129.G2) * T2_Xe129;
% T2_Xe131 = 1 / prm.WP.Xe131.G2;
% T2_Xe131_err = (prm.WP.Xe131.G2err / prm.WP.Xe131.G2) * T2_Xe131;


meas.FID_T2_129 = T2_Xe129;
meas.FID_T2_129_err = T2_Xe129_err;
meas.FID_T2_131 = T2_Xe131; 
meas.FID_T2_131_err = T2_Xe131_err;

meas.FID_B_129 = sine3fit.A129 / prm.WP.sensitivity;
meas.FID_B_129_err = fitresult_cint(1) / prm.WP.sensitivity; 
meas.FID_B_131 = sine3fit.A131 / prm.WP.sensitivity;
meas.FID_B_131_err = fitresult_cint(5) / prm.WP.sensitivity;
% 
% meas.FID_a85 = a85;
% meas.FID_a85_err = a85_err;
% meas.FID_a87 = a87;
% meas.FID_a87_err = a87_err;
% meas.FID_f85 = f85;
% meas.FID_f85_err = f85_err;
% meas.FID_f87 = f87;
% meas.FID_f87_err = f87_err;
% meas.FID_phi85 = phi85;
% meas.FID_phi85_err = phi85_err;
% meas.FID_phi87 = phi87;
% meas.FID_phi87_err = phi87_err;


% logger.log('Xenon FID measurement done')
% if ~flg.is_part_of_loop
%     disp('[+] Xenon FID measurement done');
% end


% %% Analyze
% % set local variables
% Magnetometer = v(1,:); 
% 
% logger.log('Xenon FID data processing.')
% disp('[+] Xenon FID data processing');
% 
% % estimating amplitudes for coarse fitting
% if flg.measure_magnetometer_sensitivity
% prm.WP.sensitivity = abs(2 * trapz(t, Magnetometer .* exp(-2i * pi *...
%     prm.WP.Vy_modulation.frequency * t))  / (t(end) - t(1))) /...
%     prm.WP.Vy_modulation.Bamp;
% else
%     prm.WP.sensitivity = 1;
% end
% prm.WP.Xe129.B = abs(2 * trapz(t(t < 0), Magnetometer(t < 0) .* exp(-2i *...
%     pi * prm.WP.Xe129.drive_nom.p(1) * t(t < 0))) / (max(t(t < 0)) -...
%     min(t(t < 0))) / prm.WP.sensitivity);
% prm.WP.Xe131.B = abs(2 * trapz(t(t < 0), Magnetometer(t < 0) .* exp(-2i *...
%     pi * prm.WP.Xe131.drive_nom.p(1) * t(t < 0))) / (max(t(t < 0)) -...
%     min(t(t < 0))) / prm.WP.sensitivity);
% sine3model = fittype([...
%     num2str(prm.WP.Xe129.B * prm.WP.sensitivity) '*exp(-G129*x)*cos(2*pi*f129*x+p129) + ' ...
%     num2str(prm.WP.Xe131.B * prm.WP.sensitivity) '*exp(-G131*x)*cos(2*pi*f131*x+p131) + ' ...
%     num2str(prm.WP.Vy_modulation.Bamp * prm.WP.sensitivity) '*cos(2*pi*' num2str(prm.WP.Vy_modulation.frequency) '*x+py) + dc'],...
%     'coefficients',{'G129','f129','p129','G131','f131','p131','py','dc'});
% 
% 
% % A time buffer for when the signal start to decay
% t0_for_fit = 5e-3 + 0*Tdelay;
% 
% % coarse fit using amplitudes estimation
% sine3fit = fit(t(t > t0_for_fit).', Magnetometer(t > t0_for_fit).', sine3model, 'StartPoint',...
%     [0.12, prm.WP.Xe129.drive_nom.p(1), 0.5, 0.12, prm.WP.Xe131.drive_nom.p(1), 0.5, 0.5, mean(Magnetometer)]);
% sine3model = fittype([...
%     'A129*exp(-G129*x)*cos(2*pi*f129*x+p129) + ' ...
%     'A131*exp(-G131*x)*cos(2*pi*f131*x+p131) + ' ...
%     'Ay*cos(2*pi*' num2str(prm.WP.Vy_modulation.frequency) '*x+py) + dc'],...
%     'coefficients',{'A129','G129','f129','p129','A131','G131','f131','p131','Ay','py','dc'});
% 
% % fine fit using coarse fit results
% sine3fit = fit(t(t > t0_for_fit).', Magnetometer(t > t0_for_fit).', sine3model, 'StartPoint',...
%     [prm.WP.Xe129.B * prm.WP.sensitivity, sine3fit.G129, sine3fit.f129, sine3fit.p129,...
%      prm.WP.Xe131.B * prm.WP.sensitivity, sine3fit.G131, sine3fit.f131, sine3fit.p131,...
%      prm.WP.Vy_modulation.Bamp * prm.WP.sensitivity, sine3fit.py, sine3fit.dc]);
% 
%  fitresult_int = confint(sine3fit, erf(1 / sqrt(2)));
%  fitresult_int = diff(fitresult_int,1)/2;
% 
%  
% % save relevant parameters 
% prm.WP.Xe129.G2 = sine3fit.G129; 
% prm.WP.Xe129.G2err = fitresult_int(2); 
% prm.WP.Xe129.fFID = sine3fit.f129;
% prm.WP.Xe131.G2 = sine3fit.G131; 
% prm.WP.Xe131.G2err = fitresult_int(6); 
% 
% prm.WP.Xe131.fFID = sine3fit.f131;
% 
% T2_Xe129 = 1 / prm.WP.Xe129.G2;
% T2_Xe129_err = (prm.WP.Xe129.G2err / prm.WP.Xe129.G2) * T2_Xe129;
% T2_Xe131 = 1 / prm.WP.Xe131.G2;
% T2_Xe131_err = (prm.WP.Xe131.G2err / prm.WP.Xe131.G2) * T2_Xe131;
% 
% % plot Xenon FID fit
% figure;
% hold all
% grid on;
% cmap = lines;
% t_disp = t;%t(1:round(length(t) / 16000):end);
% Magnetometer_disp = Magnetometer(1:round(length(t) / 16000):end);
% 
% dt = diff(t_disp(1:2));
% t_fit = t_disp(1):(dt/10):t_disp(end);
% plot(t, Magnetometer,'.-','color',cmap(1,:))
% plot(t_fit(t_fit > t0_for_fit), feval(sine3fit, t_fit(t_fit > t0_for_fit)), '--','color',cmap(2,:));
% % plot(sine3fit,t(t > t0_for_fit).', Magnetometer(t > t0_for_fit).');
% if flg.measure_magnetometer_sensitivity
%     magsen_str = ['magnetometer sensitivity = ' num2str(prm.WP.sensitivity / 1e3) 'V/mG\n'];
% else
%     magsen_str = 'magnetometer sensitivity = NA\n';
% end
% ttl_str = sprintf([magsen_str...
% 'B[^{129}Xe]=' num2str(prm.WP.Xe129.B * 1e6, 3) ' \\muG , T_2[^{129}Xe]=' numerr2str(T2_Xe129,T2_Xe129_err) ' s , '...
% 'B[^{131}Xe]=' num2str(prm.WP.Xe131.B * 1e6, 3) ' \\muG , T_2[^{131}Xe]=' numerr2str(T2_Xe131,T2_Xe131_err) ' s\n'...
% 'f[^{129}Xe]=' num2str(prm.WP.Xe129.drive_nom.p(1)) ' Hz , f[^{131}Xe]=' num2str(prm.WP.Xe131.drive_nom.p(1)) ' Hz']);
% title(ttl_str); 
% xlabel('time (s)'); 
% ylabel('signal (V)'); 
% xlim([-1, 9] * Tscale_Xe_FID);
% 
% logger.log('Xenon FID measurement done')
% disp('[+] Xenon FID measurement done');
% 

