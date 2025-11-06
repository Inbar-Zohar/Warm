% ========================================================================
%   This script is for measuring Xenon FID, computing the T2s and resonance 
%   frequncies from the data. The code assumes Bz, Bx and By are set, 
%   (NMRs are working) and the PIDs are disabled. 
%   The signal is measured from the scope's channel 2 which is
%   assumed to be the signal after the ESR LIA.
%
% ========================================================================
%
%   TODO:
%   1. 
% ========================================================================
% 
% if isfield(DorB,'Xe129_PID'); DorB.Xe129_PID.apply('setOutputModeMAN'); end
% if isfield(DorB,'Xe131_PID'); DorB.Xe131_PID.apply('setOutputModeMAN'); end
% 
% DorB.Bx.apply('DC',prm.WP.fields.Bx.B0);
% DorB.By.apply('DC',prm.WP.fields.By.B0);
% DorB.Xe129Drive.apply('Syncio',1);
% DorB.Xe131Drive.apply('Syncio',1);
% DorB.ESRDrive.apply('OutputON');

Tscale_XeFID = 10;
Tscale_SteadyState = 0.1;

% Init values
Vrange_XeFID = 1*[-1 1];
Vrange_magnetometer = 0.02*[-1 1];

if flg.use_DAQ_4_Xe_FID
    error('using DAQ for Xe FID is not supported yet')
else
    FID_scope = scope2;
    FID_ch = 2;
    
    % FID_scope.Reset
    % FID_scope.OperationComplete;
    % FID_scope.setChDisp(FID_ch);
    FID_scope.TrigSource('EXT');
    FID_scope.setTref('CENTer');   
    FID_scope.setTscale(Tscale_SteadyState)
    FID_scope.setTdelay(0);    
    FID_scope.setVscale(FID_ch,Vrange_XeFID)

end

init_NMR_stt = DorB.Xe129Drive.apply('Output') || DorB.Xe131Drive.apply('Output');

% measure magnetometer responsivity before FID If NMR signals are off to start with
if flg.measure_magnetometer_responsivity
    if ~init_NMR_stt
        logger.log('Measuring magnetometer responsivity');
        % apply sine modulation to By far from NMR resonance
        DorB.By.apply('Sin', 220, 5e-3, 0, prm.WP.fields.By.B0);
        [~,wvfm] = DorB.By.apply('Output',[],true);
        prm.WP.Vy_modulation.frequency = wvfm{1}.p(1);
        prm.WP.Vy_modulation.Vpp = wvfm{1}.p(2);
        prm.WP.Vy_modulation.Amplitude = prm.WP.Vy_modulation.Vpp; % Needed for legacy?
        prm.WP.Vy_modulation.Bamp = prm.WP.Vy_modulation.Vpp * prm.calib.coils.dBy_dVy / 2;
        
        % Set scope for magnetometeric measurement
        FID_scope.TrigSource('EXT');
        
        FID_scope.setVrange(FID_ch,Vrange_magnetometer)
        FID_scope.autoVrange(FID_ch,Tscale_SteadyState,[],1);
        
        FID_scope.Single; 
        FID_scope.ForceTrig;
        FID_scope.readyToRead;

        [t_mag, v_mag] = FID_scope.Read(FID_ch);
        
        DorB.By.apply('DC', prm.WP.fields.By.B0);
    end
    
else
    prm.WP.Vy_modulation.frequency = 0;
    prm.WP.Vy_modulation.Bamp = 0;
end

% Turn on NMR signals if neccesary and let the signal stabilize
logger.log('Letting Xe NMR signal stabilize');
if ~init_NMR_stt
    DorB.Xe129Drive.apply('OutputON');
    DorB.Xe131Drive.apply('OutputON');
end
pause(5 * Tscale_XeFID)

% Setup scope for Xe NMR measurement
FID_scope.autoVrange(FID_ch,Tscale_SteadyState,[],1);
% FID_scope.setVrange(FID_ch,Vrange_XeFID);
Vrange_XeFID = FID_scope.getVrange(FID_ch);

% Measure steady-state NMR
logger.log('Measuring steady-state Xe NMR signal');
FID_scope.Single;
FID_scope.ForceTrig;
FID_scope.readyToRead;
[t_NMR, v_NMR] = FID_scope.Read(FID_ch);

% Setup scope for Xe FID measurement
Tdelay = 4.0 * Tscale_XeFID;
FID_scope.setTscale(Tscale_XeFID);
FID_scope.setTdelay(Tdelay);
FID_scope.Single;

% Start Xe FID measurement
logger.log('Switching off NMRs for FID measurement');
FID_scope.ForceTrig;
pause(1 * Tscale_XeFID);
DorB.Xe129Drive.apply('OutputOFF');
DorB.Xe131Drive.apply('OutputOFF');
% pause(5 * Tscale_XeFID);
FID_scope.readyToRead;

[t_FID, v_FID] = FID_scope.Read(FID_ch);

% measure magnetometer responsivity after FID If NMR signals are on to start with
if flg.measure_magnetometer_responsivity && init_NMR_stt
    pause(1 * Tscale_XeFID);

    logger.log('Measuring magnetometer responsivity');

    DorB.By.apply('Sin', 220, 5e-3, 0, prm.WP.fields.By.B0);
    [~,wvfm] = DorB.By.apply('Output',[],true);
    prm.WP.Vy_modulation.frequency = wvfm{1}.p(1);
    prm.WP.Vy_modulation.Vpp = wvfm{1}.p(2);
    prm.WP.Vy_modulation.Bamp = prm.WP.Vy_modulation.Vpp * prm.calib.coils.dBy_dVy / 2;
    
    % Set scope for magnetometeric measurement
    FID_scope.setTscale(Tscale_SteadyState);
    FID_scope.setTdelay(0);

    FID_scope.TrigSource('EXT');
    
    FID_scope.autoVrange(FID_ch,Tscale_SteadyState,[],1);
%         FID_scope.setVrange(FID_ch,Vrange_magnetometer);

    FID_scope.Single;
    FID_scope.ForceTrig;
    FID_scope.readyToRead;
    
    [t_mag, v_mag] = FID_scope.Read(FID_ch);
    
    DorB.By.apply('DC', prm.WP.fields.By.B0);
    FID_scope.setTscale(Tscale_XeFID);
end

% Leave signals and scope in a steady state
logger.log('Switching on NMRs');
DorB.Xe129Drive.apply('OutputON');
DorB.Xe131Drive.apply('OutputON');

FID_scope.setTref('LEFT'); 
FID_scope.setTdelay(0);
FID_scope.setVrange(FID_ch,Vrange_XeFID);
FID_scope.TrigSource('LINE'); 
FID_scope.Run;



%% Save data
meas.FID_Xe_t = t_FID;
meas.FID_Xe_v = v_FID(1,:);
meas.FID_Xe_SS_t = t_NMR;
meas.FID_Xe_SS_v = v_NMR(1,:);

clear t_FID v_FID t_NMR v_NMR

if flg.measure_magnetometer_responsivity
    meas.FID_mag_t = t_mag;
    meas.FID_mag_v = v_mag(1,:);
    clear t_mag v_mag
end

if 0
    meas = create_fig(meas,10);
    subplot(311);
    hold all
    grid on;
    plot(meas.FID_mag_t,meas.FID_mag_v);
 
    subplot(312);
    hold all
    grid on;
    plot(meas.FID_Xe_SS_t,meas.FID_Xe_SS_v);

    subplot(313);
    hold all
    grid on;
    plot(meas.FID_Xe_t,meas.FID_Xe_v);
end


logger.log('Xenon FID data processing.')

%% analyze results

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

if flg.measure_magnetometer_responsivity
    % Get magnetometer responsivity from the measurement with 
    % (almost?) no NMR signal
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
    
else
    Amag = 0;
    mag_sens = 1;
    mag_sens_err = 0;
end

%%
FID_Xe_t = meas.FID_Xe_t;
FID_Xe_v = meas.FID_Xe_v;

% Get Xe magnetization estimate from steady-state NMR signal

% A time buffer for steady state
t0_for_fit_NMR = -25e-3;
FID_Xe_NMR_t = meas.FID_Xe_t(meas.FID_Xe_t < t0_for_fit_NMR);
FID_Xe_NMR_v = meas.FID_Xe_v(meas.FID_Xe_t < t0_for_fit_NMR);

dc_guess = mean(FID_Xe_NMR_v);
A_guess = (max(FID_Xe_NMR_v) - min(FID_Xe_NMR_v))/2;

NMR_fittype = fittype(@(A129,A131,p129,p131,pmag,dc,x)...
    sine3_ss_func(A129,A131,0*Amag,f_129_drive,f_131_drive,f_mag,p129,p131,pmag,dc,x));
NMR_fit = fit(FID_Xe_NMR_t(:), FID_Xe_NMR_v(:), NMR_fittype, ...
    'StartPoint', [A_guess/2 A_guess/2 0 0 0 dc_guess],...
    'Lower',[0 0 -pi -pi -pi -Inf], ...
    'Upper', [Inf Inf 2*pi 2*pi 2*pi Inf]);

% Finally do FID analysis

% A time buffer for when the signal start to decay
t0_for_fit_FID = 25e-3;
FID_Xe_t_for_fit = FID_Xe_t(FID_Xe_t > t0_for_fit_FID); 
FID_Xe_v_for_fit = FID_Xe_v(FID_Xe_t > t0_for_fit_FID); 
max_t = FID_Xe_t_for_fit(end);

% Try to find get precession frequencies from 
[~, ~, ~, fft_freqs] = fft1(FID_Xe_t_for_fit(:), FID_Xe_v_for_fit(:),1,0);
if sum(abs(fft_freqs - f_129_drive) < 1) == 1
    f_129_guess = fft_freqs(abs(fft_freqs - f_129_drive) < 1);
else
    f_129_guess = f_129_drive;
end
if sum(abs(fft_freqs - f_131_drive) < 1) == 1
    f_131_guess = fft_freqs(abs(fft_freqs - f_131_drive) < 1);
else
    f_131_guess = f_131_drive;
end
% f_129_guess = f_129_drive;
% f_131_guess = f_131_drive;


% coarse fit using amplitudes estimation from steady-state NMR fit
% FID_fittype_coarse = fittype(@(f129,f131,p129,p131,pmag,G129,G131,x)...
%     sine3_decay_func(NMR_fit.A129,NMR_fit.A131,0,...
%     f129,f131,f_mag,p129,p131,0,G129,G131,NMR_fit.dc,0,x));
% FID_fit_coarse = fit(FID_Xe_t_for_fit(:), FID_Xe_v_for_fit(:), FID_fittype_coarse, ...
%     'StartPoint', [f_129_drive f_131_drive 0 0 2/max_t 2/max_t]);%,...
% %     'Lower',[0 0 0 0 0 0], ...
% %     'Upper', [Inf Inf 2*pi 2*pi Inf Inf]);

FID_fittype_coarse_no_freq = fittype(@(p129,p131,pmag,G129,G131,x)...
    sine3_decay_func(NMR_fit.A129,NMR_fit.A131,0,...
    f_129_guess,f_131_guess,f_mag,p129,p131,pmag,G129,G131,NMR_fit.dc,0,x));
FID_fit_coarse_no_freq = fit(FID_Xe_t_for_fit(:), FID_Xe_v_for_fit(:), FID_fittype_coarse_no_freq, ...
    'StartPoint', [0 0 0 2/max_t 2/max_t],...
    'Lower',[-pi -pi -pi 0 0], ...
    'Upper', [2*pi 2*pi 2*pi Inf Inf]);

% fine fit using coarse fit results
FID_fittype = fittype(@(A129,A131,Amag,f129,f131,p129,p131,pmag,G129,G131,dc,x)...
    sine3_decay_func(A129,A131,Amag,...
    f129,f131,f_mag,p129,p131,pmag,G129,G131,dc,0,x));
FID_fit = fit(FID_Xe_t_for_fit(:), FID_Xe_v_for_fit(:), FID_fittype, ...
    'StartPoint', [NMR_fit.A129,NMR_fit.A129,0,...
    f_129_guess f_131_guess ...
    FID_fit_coarse_no_freq.p129 FID_fit_coarse_no_freq.p131 0 ...
    FID_fit_coarse_no_freq.G129 FID_fit_coarse_no_freq.G131 NMR_fit.dc],...
    'Lower',[0 0 0 0 0 -pi -pi -pi 0 0 -Inf], ...
    'Upper', [Inf Inf Inf Inf Inf 2*pi 2*pi 2*pi Inf Inf Inf ]);

fitresult_cint = confint(FID_fit, erf(1 / sqrt(2)));
fitresult_cint = diff(fitresult_cint,1)/2;

% Extract interesting params from the fit results
G129 = FID_fit.G129;
G129_err = fitresult_cint(9);
G131 = FID_fit.G131;
G131_err = fitresult_cint(10);

T2_Xe129 = 1 / G129;
T2_Xe129_err = (G129_err / G129) * T2_Xe129;
T2_Xe131 = 1 / G131;
T2_Xe131_err = (G131_err / G131) * T2_Xe131;

A129 = FID_fit.A129;
A129_err = fitresult_cint(1);
A131 = FID_fit.A131;
A131_err = fitresult_cint(2);

B129 = A129/mag_sens;
B129_err = sqrt((A129_err/mag_sens)^2 + (B129*mag_sens_err/mag_sens)^2);
B131 = A131/mag_sens;
B131_err = sqrt((A131_err/mag_sens)^2 + (B131*mag_sens_err/mag_sens)^2);

fFID_129 = FID_fit.f129;
fFID_129_err = fitresult_cint(5);
fFID_131 = FID_fit.f131;
fFID_131_err = fitresult_cint(6);


%% save relevant parameters 
meas.FID_T2_129 = T2_Xe129;
meas.FID_T2_129_err = T2_Xe129_err;
meas.FID_T2_131 = T2_Xe131; 
meas.FID_T2_131_err = T2_Xe131_err;

meas.FID_B_129 = B129;
meas.FID_B_129_err = B129_err; 
meas.FID_B_131 = B131;
meas.FID_B_131_err = B131_err;

meas.FID_f_129 = fFID_129;
meas.FID_f_129_err = fFID_129_err; 
meas.FID_f_131 = fFID_131;
meas.FID_f_131_err = fFID_131_err;

meas.FID_mag_sens = mag_sens;
meas.FID_mag_sens_err = mag_sens_err;

prm.WP.responsivity = mag_sens;

%% plot results 

% plot Xenon FID fit
cmap = lines;


typical_data_dt = 0.25e-3;
fit_to_data_factor = 2;

if flg.measure_magnetometer_responsivity
    %%
    meas = create_fig(meas,6);
    clf;
    set(gcf,'position', [300 300 900 300]);
    
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

    magsen_str = ['magnetometer responsivity = ' ...
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
    
else
    magsen_str = '';
end

downsample_factor_data = ceil(typical_data_dt/diff(FID_Xe_t(1:2)));
inds_disp_data = 1:downsample_factor_data:length(FID_Xe_t);
t_disp_data = FID_Xe_t(inds_disp_data);
% v_disp_data = FID_Xe_v(inds_disp_data);
v_disp_data = decimate(FID_Xe_v,downsample_factor_data);
t_disp_data = t_disp_data(1:length(v_disp_data));
t_disp_fit = linspace(t_disp_data(1),t_disp_data(end),length(t_disp_data)*5);

t_disp_fit_NMR =  t_disp_fit(t_disp_fit < t0_for_fit_NMR);
fit_disp_NMR = feval(NMR_fit, t_disp_fit_NMR);

t_disp_fit_FID = t_disp_fit(t_disp_fit > t0_for_fit_FID);
fit_disp_FID = feval(FID_fit, t_disp_fit_FID);
fit_envelope_129 = FID_fit.A129 *exp(-FID_fit.G129*t_disp_fit_FID);
fit_envelope_131 = FID_fit.A131 *exp(-FID_fit.G131*t_disp_fit_FID);
fit_envelope_fixed = FID_fit.dc + [1;-1]*(FID_fit.Amag + 0*t_disp_fit_FID) ;

fit_envelope_129_disp = fit_envelope_fixed + [1;-1]*fit_envelope_129;
fit_envelope_131_disp = fit_envelope_fixed + [1;-1]*fit_envelope_131;
fit_envelope_tot_disp = fit_envelope_fixed + [1;-1]*(fit_envelope_129+fit_envelope_131);

% Display envelope at magneteometer drive frequency only if it's large enough
if abs(FID_fit.Amag) < 1e-4; fit_envelope_fixed = fit_envelope_fixed*nan; end

yl = FID_fit.dc + 1.1 * [-1 1]*(FID_fit.Amag + FID_fit.A129 +  FID_fit.A131);
yl(1) = min(yl(1), min(FID_Xe_v) - 0.1 * (max(FID_Xe_v) - min(FID_Xe_v)));
yl(2) = max(yl(2), max(FID_Xe_v) + 0.1 * (max(FID_Xe_v) - min(FID_Xe_v)));

zoom_t_range = round(1.5 / fFID_131,2);
zoom_t_start = [round(mean([t_disp_data(1),t0_for_fit_NMR]),1) t0_for_fit_FID  round(FID_Xe_t(end)/3) round(0.9*FID_Xe_t(end))];
zoom_inds_data = [];
zoom_inds_fit = [];

for ind=1:length(zoom_t_start)
    this_inds = find(t_disp_data >=  zoom_t_start(ind) & t_disp_data <=  zoom_t_start(ind) + zoom_t_range);
    zoom_inds_data(ind,:) = this_inds([1 end]);
    if zoom_t_start(ind) >= t0_for_fit_FID
        this_inds = find(t_disp_fit_FID >=  zoom_t_start(ind) & t_disp_fit_FID <=  zoom_t_start(ind) + zoom_t_range);
    else
        this_inds = find(t_disp_fit_NMR >=  zoom_t_start(ind) & t_disp_fit_NMR <=  zoom_t_start(ind) + zoom_t_range);
    end
    zoom_inds_fit(ind,:) = this_inds([1 end]);
end

meas = create_fig(meas,5);
set(gcf,'position', [300 300 900 650]);

subplot(2,length(zoom_t_start),1:length(zoom_t_start))
hold all
grid on;
box on;
pos = get(gca,'position');
set(gca,'position', pos + [0 -0.09 0 0.02]);

for ind=1:length(zoom_t_start)
    patch('Xdata',zoom_t_start(ind)+zoom_t_range*[0 0 1 1], 'Ydata',yl([1 2 2 1]), ...
        'facecolor', 0.5*ones(1,3), 'facealpha', 0.3, 'edgecolor', 0.5*ones(1,3)) 

end

plot(t_disp_data, v_disp_data,'.','color',cmap(1,:))
plot(t_disp_fit_FID, fit_disp_FID, '-','color',cmap(2,:));
plot(t_disp_fit_NMR, fit_disp_NMR, '-','color',cmap(2,:));

plot(t_disp_fit_FID, fit_envelope_fixed, ':','color','k','linewidth', 2);
plot(t_disp_fit_FID, fit_envelope_129_disp, ':','color',cmap(3,:),'linewidth', 2);
plot(t_disp_fit_FID, fit_envelope_131_disp, ':','color',cmap(4,:),'linewidth', 2);
plot(t_disp_fit_FID, fit_envelope_tot_disp, '--','color','k','linewidth', 2);



ttl_str = sprintf([char(meas.DateTime_var), '\n', magsen_str,...
'B[^{129}Xe]=' numerr2str(B129*1e6,B129_err*1e6) ' \\muG , T_2[^{129}Xe]=' numerr2str(T2_Xe129,T2_Xe129_err) ' s , '...
'B[^{131}Xe]=' numerr2str(B131*1e6,B131_err*1e6) ' \\muG , T_2[^{131}Xe]=' numerr2str(T2_Xe131,T2_Xe131_err) ' s\n'...
'f[^{129}Xe]=' numerr2str(fFID_129,fFID_129_err) ' Hz (drive at ' num2str(f_129_drive) 'Hz), '...
'f[^{131}Xe]=' numerr2str(fFID_131,fFID_131_err) ' Hz (drive at ' num2str(f_131_drive) 'Hz)']);
title(ttl_str); 
xlabel('time (s)'); 
ylabel('signal (V)'); 
xlim(FID_Xe_t([1 end]));
ylim(yl);

yyaxis right;
xlim(FID_Xe_t([1 end]));
set(gca,'YColor','black');
ylabel('signal (\muG)'); 
ylim((yl-FID_fit.dc) /mag_sens *1e6);

yyaxis left;

for ind=1:length(zoom_t_start)
    subplot(2,length(zoom_t_start),length(zoom_t_start)+ind)
    hold all
    grid on;
    box on;
    pos = get(gca,'position');
    set(gca,'position', pos - [0 0.04 0 0]);

    plot((t_disp_data(zoom_inds_data(ind,1):zoom_inds_data(ind,2))-zoom_t_start(ind))*1e3, ...
        v_disp_data(zoom_inds_data(ind,1):zoom_inds_data(ind,2)),...
        '.','color',cmap(1,:),'markersize',8)
    
    if zoom_t_start(ind) >= t0_for_fit_FID
        plot((t_disp_fit_FID(zoom_inds_fit(ind,1):zoom_inds_fit(ind,2))-zoom_t_start(ind))*1e3, ...
            fit_disp_FID(zoom_inds_fit(ind,1):zoom_inds_fit(ind,2)), ...
            '-','color',cmap(2,:));
        plot((t_disp_fit_FID(zoom_inds_fit(ind,1):zoom_inds_fit(ind,2))-zoom_t_start(ind))*1e3, ...
            fit_envelope_fixed(:,zoom_inds_fit(ind,1):zoom_inds_fit(ind,2)), ...
            ':','color','k','linewidth', 2);
        plot((t_disp_fit_FID(zoom_inds_fit(ind,1):zoom_inds_fit(ind,2))-zoom_t_start(ind))*1e3, ...
            fit_envelope_129_disp(:,zoom_inds_fit(ind,1):zoom_inds_fit(ind,2)), ...
            ':','color',cmap(3,:),'linewidth', 2);
        plot((t_disp_fit_FID(zoom_inds_fit(ind,1):zoom_inds_fit(ind,2))-zoom_t_start(ind))*1e3, ...
            fit_envelope_131_disp(:,zoom_inds_fit(ind,1):zoom_inds_fit(ind,2)), ...
            ':','color',cmap(4,:),'linewidth', 2);
        plot((t_disp_fit_FID(zoom_inds_fit(ind,1):zoom_inds_fit(ind,2))-zoom_t_start(ind))*1e3, ...
            fit_envelope_tot_disp(:,zoom_inds_fit(ind,1):zoom_inds_fit(ind,2)), ...
            '--','color','k','linewidth', 2);
    else
        plot((t_disp_fit_NMR(zoom_inds_fit(ind,1):zoom_inds_fit(ind,2))-zoom_t_start(ind))*1e3, ...
            fit_disp_NMR(zoom_inds_fit(ind,1):zoom_inds_fit(ind,2)), ...
            '-','color',cmap(2,:));
    end
    xlim([0, zoom_t_range*1e3]);
    ylim(yl);

    if abs(zoom_t_start(ind)) >= 0.1
        t0str = [num2str(zoom_t_start(ind)) 's'];
    else
        t0str = [num2str(zoom_t_start(ind)*1e3) 'ms'];
    end
    if zoom_t_start(ind) < 0
        t0str = ['(' t0str ')'];
    end
    xlabel(['t-' t0str ' (ms)']);
        
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
        ylim((yl-FID_fit.dc) /mag_sens *1e6);
        
    end
   
end
%%

