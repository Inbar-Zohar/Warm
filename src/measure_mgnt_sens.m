%measure magnetometer sens via WR with shot drivers.

%% Params
LIA_gain = 10/500e-3;
Tscale_XeFID = 5;
Tscale_SteadyState = 5e-3;
By_amppp = 50e-3;

% Init values
Vrange_XeFID = 1*[-1 1];
Vrange_magnetometer = 0.02*[-1 1];


% Turn off Xe drives
init_NMR_stt = DorB.Xe129Drive.apply('Output') || DorB.Xe131Drive.apply('Output');

if init_NMR_stt %turning off NMR signal if needed
    logger.log('Switching off Xenon drives');
    DorB.Xe129Drive.apply('OutputOFF');
    DorB.Xe131Drive.apply('OutputOFF');
    logger.log('Waiting for Xenon to decay');
    pause(5 * Tscale_XeFID)
end

%% Measure magnetometer sensitivity
if prm.mgnt_meas_magsens
    logger.log('Measuring magnetometer sensitivity');
    
    % apply sine modulation to By far from NMR resonance
    DorB.By.apply('Sin', 220, By_amppp, 0, prm.WP.fields.By.B0);
    [~,wvfm] = DorB.By.apply('Output',[],true);
    prm.WP.Vy_modulation.frequency = wvfm{1}.p(1);
    prm.WP.Vy_modulation.Vpp = wvfm{1}.p(2);
    prm.WP.Vy_modulation.Amplitude = prm.WP.Vy_modulation.Vpp;
    prm.WP.Vy_modulation.Bamp = prm.WP.Vy_modulation.Vpp * prm.calib.coils.dBy_dVy / 2;

    % Set scope for magnetometeric measurement
    FID_scope = scope2;
    FID_scope.setTscale(Tscale_SteadyState);
    FID_scope.setTdelay(0);
    
    FID_ch = 2;

    FID_scope.TrigSource('EXT');

    %FID_scope.setVrange(FID_ch,Vrange_magnetometer)
    FID_scope.setVrange(FID_ch,5*[-1 1])
    FID_scope.autoVrange(FID_ch,Tscale_SteadyState,[],1);

    FID_scope.Single; 
    FID_scope.ForceTrig;
    FID_scope.readyToRead;

    [t_mag, v_mag] = FID_scope.Read(FID_ch);

    DorB.By.apply('DC', prm.WP.fields.By.B0);

    pause(5 * Tscale_SteadyState)

    %%saving data

     meas.FID_mag_t = t_mag;
     meas.FID_mag_v = v_mag(1,:);
     clear t_mag v_mag


     % Common definitions
    sine3_ss_func = @(A129,A131,Amag,f129,f131,fmag,p129,p131,pmag,dc,x) ...
        A129.*cos(2*pi*f129*x+p129) + ...
        A131.*cos(2*pi*f131*x+p131) + ...
        Amag.*cos(2*pi*fmag*x+pmag) + dc;

    sine3_decay_func = @(A129,A131,Amag,f129,f131,fmag,p129,p131,pmag,G129,G131,dc,x0,x) ...
        A129.*cos(2*pi*f129*x+p129).*exp(-G129*(x-x0)) + ...
        A131.*cos(2*pi*f131*x+p131).*exp(-G131*(x-x0)) + ...
        Amag.*cos(2*pi*fmag*x+pmag) + dc;


f_129_drive = prm.WP.Xe129.freq;
f_131_drive = prm.WP.Xe131.freq;


     FID_mag_t = meas.FID_mag_t;
     FID_mag_v = meas.FID_mag_v;

    f_mag = prm.WP.Vy_modulation.frequency;
    dc_guess = mean(FID_mag_v);
    A_guess = (max(FID_mag_v) - min(FID_mag_v))/2;

    % Try to fit assuming all 3 frequencies exist
    mag_fittype = fittype(@(A129,A131,Amag,f129,f131,p129,p131,pmag,dc,x)...
        sine3_ss_func(A129,A131,Amag,f129,f131,f_mag,p129,p131,pmag,dc,x));
    mag_fit = fit(FID_mag_t(:), FID_mag_v(:), mag_fittype, ...
        'StartPoint', [0 0 A_guess f_129_drive f_131_drive 0 0 0 dc_guess],...
        'Lower',[0 0 0 0 0 -pi -pi -pi -Inf], ...
        'Upper', [Inf Inf Inf Inf Inf 2*pi 2*pi 2*pi Inf]);

    fitresult_cint = confint(mag_fit, erf(1 / sqrt(2)));
    fitresult_cint = diff(fitresult_cint,1)/2;
    Amag_err = fitresult_cint(3);

    % If DC uncertainty too large, it means the fit is underconstrained. We
    % can probably fix it by forcing the 129 contribution to be zero.
    if fitresult_cint(end) > A_guess
        mag_fittype = fittype(@(A131,Amag,f131,p131,pmag,dc,x)...
            sine3_ss_func(0,A131,Amag,0,f131,f_mag,p131,p131,pmag,dc,x));        
        mag_fit = fit(FID_mag_t(:), FID_mag_v(:), mag_fittype, ...
            'StartPoint', [0 A_guess f_131_drive 0 0 dc_guess],...
            'Lower',[0 0 0 -pi -pi -Inf], ...
            'Upper', [Inf Inf Inf 2*pi 2*pi Inf]);
        fitresult_cint = confint(mag_fit, erf(1 / sqrt(2)));
        fitresult_cint = diff(fitresult_cint,1)/2;
        Amag_err = fitresult_cint(2); 
    end

    Amag = mag_fit.Amag;    
    dc_mag = mag_fit.dc;
    mag_sens = mag_fit.Amag / prm.WP.Vy_modulation.Bamp;
    mag_sens_err = Amag_err / prm.WP.Vy_modulation.Bamp;

    meas.FID_mag_sens = mag_sens;
    meas.FID_mag_sens = mag_sens_err;

    prm.WP.sensitivity = mag_sens;
end

%% measuring ESR LIA noise floor
logger.log('Measuring ESR noise floor');
% Initial setup
prm.to_plot = prm.mgnt_points_to_plot;
prm.readtimeS = floor(prm.mgnt_exp_time_s);

prm.cDAQ_channels = [2]; %Only ESR

% get values from instruments
[~, sys_data.DATATESR] = AG3.Output(1, [], true);
sys_data.DATATESR = sys_data.DATATESR{1}.p(1:2);
prm.fresESR = sys_data.DATATESR(1);  %[Hz]
prm.AmpESR = [sys_data.DATATESR(2)]; %[Vpp]

if prm.mgnt_Fs >51200
    prm.mgnt_Fs = 51200;
    warning('Fs was set to its maximal value of 51200 samp/sec');
end

% read DAQ
pause(prm.PAUSE_TIME_SEC);
cDAQ.setFs(prm.mgnt_Fs);
cDAQ.setDuration(prm.readtimeS + prm.PAUSE_TIME_SEC);
% disp('[+] reading DAQ...');
[tpr, vpr] = cDAQ.measureForeground;
logger.log('[+] done reading DAQ');

% Turn on Xe back
if turn_on_Xe
DorB.Xe129Drive.apply('OutputON');
DorB.Xe131Drive.apply('OutputON');
end

% the data
if numel(vpr)*8>7.078051776e+09/10*2 % larger than 2 GB
    vpr_tmp = vpr(1:10:end,:);
    vpr_sngl = single(vpr_tmp);
    meas.tp_wrld_rot = tpr(1:10:end)';
    meas.vp_wrld_rot = vpr_sngl';
    clear vpr_tmp vpr_sngl
else
    meas.tp_wrld_rot = tpr';
    meas.vp_wrld_rot = vpr';
end
clear sys_data

ESR_LIA_out   = meas.vp_wrld_rot(1, :);
meas.ESR_LIA_out = ESR_LIA_out(1:end)';
meas.time_s = meas.tp_wrld_rot(1:end);

logger.log('Measurents done');

%% plot data
f_cut_off = [5e1 0.8e3];


%mag sens
if (prm.mgnt_meas_magsens)
    cmap = lines;
    meas = create_fig(meas,18);
    %myfig(3)
    clf;
    % set(gcf,'position', [300 300 900 300]);

    typical_data_dt = 0.25e-3;
    downsample_factor_data = ceil(typical_data_dt/diff(FID_mag_t(1:2)));
    inds_disp_data = 1:downsample_factor_data:length(FID_mag_t);
    t_disp_data = FID_mag_t(inds_disp_data);
    v_disp_data = FID_mag_v(inds_disp_data);

    t_disp_fit = linspace(t_disp_data(1),t_disp_data(end),length(t_disp_data)*5);
    fit_disp_mag = feval(mag_fit, t_disp_fit);
    zoom_t_range = 20e-3;

    yl = mag_fit.dc + 1.1 * [-1 1]*(mag_fit.Amag + mag_fit.A131);
    yl(1) = min(yl(1), min(v_disp_data) - 0.1 * (max(v_disp_data) - min(v_disp_data)));
    yl(2) = max(yl(2), max(v_disp_data) + 0.1 * (max(v_disp_data) - min(v_disp_data)));

    sp1 = subplot(1,3,[1 2]); cla;
    hold all
    grid on;
    box on;
    patch('Xdata',zoom_t_range*[0 0 1 1], 'Ydata',yl([1 2 2 1]), ...
        'facecolor', 0.5*ones(1,3), 'facealpha', 0.3, 'edgecolor', 0.5*ones(1,3)) 
    plot(t_disp_data, v_disp_data,'.','color',cmap(1,:))
    plot(t_disp_fit, fit_disp_mag, '-','color',cmap(2,:));
    xlabel('time (s)');
    ylabel('signal (V)'); 
    ylim(yl)

    magsen_str = ['magnetometer sensitivity = ' ...
        numerr2str(mag_sens/1e3,mag_sens_err/1e3) ' V/mG\n'];


    title(sprintf('%s, %s',char(meas.DateTime_var),magsen_str(1:end-2)));

    sp2=subplot(1,3,3);cla;
    hold all
    grid on;
    box on;

    plot(t_disp_data*1e3, v_disp_data,'.','color',cmap(1,:))
    plot(t_disp_fit*1e3, fit_disp_mag, '-','color',cmap(2,:));
    xlim([0 zoom_t_range]*1e3)
    ylim(yl)
    xlabel('time (ms)'); 
    pause(0.1)
    sp1.Position([2 4]) = sp2.Position([2 4]);
    sp1.XLabel.Position(2) = sp2.XLabel.Position(2)+1e-4; 
end

%noise floor
T = meas.time_s;
fs = 1 / mean(diff(T));
[f,P1,pks,locks,Y] = fft1(T,ESR_LIA_out,1);
[psd,fff]=pwelch(ESR_LIA_out(2:end) - mean(ESR_LIA_out(2:end)), [], [], [], fs);

ind_end = find(fff>f_cut_off(2),1);
ind_srt = find(fff>f_cut_off(1),1);


if(prm.mgnt_meas_magsens)
    srPSD_fT = sqrt(psd).*(1e-4*1e15/mag_sens);
    %find baseline [fT/srHz]
    YBfT = msbackadj(fff(ind_srt:ind_end),srPSD_fT(ind_srt:ind_end),'EstimationMethod','em'); %set windowsize and stepsize to optimom
    bgfT = srPSD_fT(ind_srt:ind_end) - YBfT;
    LIA_mgnt_noise_floor_fT_srHz_err = std(bgfT)
    LIA_mgnt_noise_floor_fT_srHz = mean(bgfT)
end

%find baseline [V/srHz]
srPSD = sqrt(psd);
YB_V = msbackadj(fff(ind_srt:ind_end),srPSD(ind_srt:ind_end),'EstimationMethod','em'); %set windowsize and stepsize to optimom
bg_V = srPSD(ind_srt:ind_end) - YB_V;
LIA_mgnt_noise_floor_V_srHz = mean(bg_V)
LIA_mgnt_noise_floor_V_srHz_err = std(bg_V)
LIA_mgnt_noise_floor_V_srHz/LIA_gain;

Nflr_V_str = ['Noise floor = ' numerr2str(LIA_mgnt_noise_floor_V_srHz*1e6,LIA_mgnt_noise_floor_V_srHz_err*1e6) ' uV/srHz'];
%plot 1 Voltage
meas = create_fig(meas,20);
%myfig(1);
clf;
set(gcf,'position', [300 300 900 300]);

%time domain
spl2_1 = subplot(3,1,1); cla;
grid on;
box on;
plot(T, ESR_LIA_out);
xlabel('Time (sec)');
ylabel('Signal ESR (V)');
title(sprintf('%s, %s',char(meas.DateTime_var),Nflr_V_str(1:end-2)));

%fft
spl2_2 = subplot(3,1,2); cla;
plot(f,P1)
xlabel('Frequency (HZ)');
ylabel('FFT ESR (V)')
box on;
grid on;

%psd
spl2_3 = subplot(3,1,3); cla;
hold on
loglog(fff,srPSD)
loglog([fff(ind_srt) fff(ind_end)],[LIA_mgnt_noise_floor_V_srHz LIA_mgnt_noise_floor_V_srHz])
loglog(fff(ind_srt:ind_end),bg_V)
hold off
xlabel('Frequency (HZ)');
ylabel('{\surd}PSD of ESR (V/{\surd}Hz)')
grid on;
box on;
set(gca,'YScale','log')
set(gca,'XScale','log')

%plot 2 fT/srHz
if (prm.mgnt_meas_magsens)
    Nflr_fT_str = ['Noise floor = ' numerr2str(LIA_mgnt_noise_floor_fT_srHz, LIA_mgnt_noise_floor_fT_srHz_err) ' \fT/srHz\n'];
    meas = create_fig(meas,30);
    %myfig(2)
    clf;
    set(gcf,'position', [300 300 900 300]);

    %time domain
    spl3_1 = subplot(3,1,1); cla;
    grid on;
    box on;
    plot(T, ESR_LIA_out.*(1e-4*1e15/mag_sens)); 
    title(sprintf('%s, %s',char(meas.DateTime_var),Nflr_fT_str(1:end-2)));
    xlabel('Time (sec)');
    ylabel('Signal ESR (fT)');

    %fft
    spl3_2 = subplot(3,1,2); cla;
    grid on;
    box on;
    plot(f,P1.*(1e-4*1e15/mag_sens))
    xlabel('Frequency (HZ)');
    ylabel('FFT ESR (fT)')

    %psd
    spl3_3 = subplot(3,1,3); cla;

    loglog(fff,srPSD_fT)
    hold on
    loglog([fff(ind_srt) fff((ind_end))],[LIA_mgnt_noise_floor_fT_srHz LIA_mgnt_noise_floor_fT_srHz])
    loglog(fff(ind_srt:ind_end),bgfT)
    hold off
    xlabel('Frequency (HZ)');
    ylabel('{\surd}PSD of ESR (fT/{\surd}Hz)')
    grid on
        box on;
    set(gca,'YScale','log')
set(gca,'XScale','log')

end

%Reset the cDAQ
get_daq;
