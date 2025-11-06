% ========================================================================
%   This script is measuring Alkali FID and computing the T2 and resonance 
%   frequncies from the data. The code assumes Bz is set, Bx and By are 
%   zeroed, BPD is balanced by a lambda plate. 
%   The signal from BPD is connected ch1 scope1. 
% ========================================================================
%
%   TODO:
%   # compute the larmour frequency from known values instead of reading it
%   from the ESR function generator
%   # set scope to averaging so that it trnafers 12bit data (check it)
%
%   Set scope scale dynamically - Roy 23/09/2021
%   Corrected the fitting routine - CF 29/09/21
%
% ========================================================================


logger.log('Alkali FID measurement initialized.')
disp('[+] Starting an Alkali FID measurement');
% Disable PID 129 131 ESR

% Disable NMR 129 and 131 and ESR
AG1.OutputOFF(1);   AG2.OutputOFF(1);   AG3.OutputOFF(1);
% Disable 

% Set  AG6 ch2 to trigger all with a square pulse
AG6.Square(2, 100, 1.5, 0, 0);  AG6.OutputON(2);

% read larmour frequency of Alkali from function generator 
[stt, wvfm] = AG3.Output(1, [], true);
prm.fresESR = wvfm{1}.p(1);

% % compute the larmour frequency from known values instead of reading the
% % ESR function generator
% prm.fresESR = prm.WP.coils.dBz_dIz * prm.WP.fields.Iz * c.gRb85;

% Set AG2 ch2 to By pulses, with ext trig
dt = 1 / 10 / prm.fresESR; % dt is shorter 10 times than larmours cycle. (NMR limit)
AG2.BurstPulse(2, 1e3, 5, dt, 0, 0, 0)

% Set scope 1 ch1 with trigger on ch1
% Configure scope to take data
chnl = 1;
Tscale = 50e-6; % sec per dev
Vscale = 0.04; % volt per dev

% set scope
curr_inst = scope1;
curr_inst.Stop;

% Coupling
curr_inst.setChCoupling(chnl, 'DC');
curr_inst.setChImpedance(chnl, 'M');

% Acquire
% curr_inst.HighRes;
% curr_inst.AverageAQN(2); % for 12bit data (not working as is)

% horizontal
curr_inst.setTscale(Tscale);
curr_inst.setTdelay(5.00 * Tscale);

% vertical offset
curr_inst.setVoffset(chnl, 0);

% setting the scope vscale dynamically
try_counter = 0;
in_once_flag = 0;
out_once_flag = 0;

while try_counter < 30 && ~in_once_flag && Vscale >= 0.002
    % setting signal to scope bounds by enlarging the scale until signal is
    % out of scope and than setting it to smaller scale until signal is in
    % again.
    
    % vertical scale
    curr_inst.setVscale(chnl, Vscale);
    curr_inst.Run;

    % Trigger
    curr_inst.TrigMode('NORM');
    curr_inst.TrigSource('EXT');
    curr_inst.TrigSlope('POS');
    curr_inst.TrigThresh(1);
    curr_inst.Single;

    % read data from scope
    curr_inst.readyToRead(10);
    [t, v, ~] = curr_inst.Read(chnl);

    in_flag = curr_inst.sigInbounds(v, chnl);
    if ~in_flag
        % signal not contained in scope, enlarge scale
        Vscale = Vscale * 1.25;
        try_counter = try_counter + 1;
        out_once_flag = 1;
        
    elseif in_flag && ~out_once_flag
        % signal is contained in scope from the start, set smaller scale
        Vscale = Vscale * 0.5;
        try_counter = try_counter + 1;
        
    elseif in_flag && out_once_flag
        in_once_flag = 1;
    end
    

end


% Release scope from Single mode
curr_inst.Run;

meas.FID_Rb_t = t;
meas.FID_Rb_v = v;


% % uncomment to plot waveform 
% figure; plot(meas.FID_Rb_t / Tscale, meas.FID_Rb_v); xlabel('time'); ylabel('amplitude') 

meas.raw_data{1,FIDindex}= t;
meas.raw_data{2,FIDindex}=v;
[~,wvfm1] = AG4.Output(1,[],true); 
[~,wvfm2] = AG4.Output(2,[],true); 
meas.raw_data{3,FIDindex}=['probe=' num2str(wvfm2{1}.p(end)) '_pump=' num2str(wvfm1{1}.p(end))];
FIDindex=FIDindex+1;
clear t v 

%% FID rubidium Analyze


logger.log('Analizing measurement results')
% % uncomment to plot waveform 
% figure; plot(meas.FID_Rb_t / Tscale, meas.FID_Rb_v); xlabel('time'); ylabel('amplitude [volts]'); hold on

% find measured signal start index
[~, ind1] = max(meas.FID_Rb_v);

% % plot truncated signal - uncomment for plotting
% figure; plot(meas.FID_Rb_t(ind1:end), meas.FID_Rb_v(ind1:end)); hold on

% set variables for fitting process
t4a = meas.FID_Rb_t(ind1:end);
v4a = meas.FID_Rb_v(ind1:end);

% Find peaks in FFT
[f, P1] = fft1(t4a, v4a);
[pks, locs, w, p] =  findpeaks(P1, f, 'Annotate', 'extents',...
    'MinPeakProminence', 0.1 * max(P1(:)), 'SortStr', 'descend');


% get envelope - T2
np2 = 1000;
[up2, lo2] = envelope(v4a, np2, 'peak');

% estimate T2 by finding the time where the envelope decreases to half
ind2 = find(up2 - up2(end) < (up2(1) - up2(end)) / 2, 1);
T2_init = t4a(ind2);
g_init = 1 / T2_init;
n4fit = 10; % the length of data for fit in units of T2
ind3 = min([numel(t4a), ind2 * n4fit]); % make sure not to overflow
[xData, yData] = prepareCurveData(t4a(1:ind3), v4a(1:ind3));

% Set up fittype and options.
fit_fnc = 'a85*sin(2*pi*f85*x+phi85)*exp(-g85*x)+a87*sin(2*pi*f87*x+phi87)*exp(-g87*x)+c';
ft = fittype(fit_fnc, 'independent', 'x', 'dependent', 'y');
opts = fitoptions('Method', 'NonlinearLeastSquares');
a_init = (up2(1) - yData(end));
opts.Display = 'Off'; %(a85,   a87,    c,                     f85,         f87,          g85,        g87,            phi85,   phi87,x)
opts.StartPoint =      [a_init a_init  yData(end)             locs(1)      locs(2)       g_init      g_init     pi        pi   ];
opts.Lower =       0.1*[a_init a_init  -100*abs(yData(end))   9*locs(1)    9*locs(2)     g_init      g_init     0         0    ];
opts.Upper =        10*[a_init a_init  abs(yData(end))        0.11*locs(1) 0.11*locs(2)  g_init      g_init  0.5*pi    0.5*pi  ];

% Fit model to data.
[fitresult, gof] = fit(xData, yData, ft, opts);

% set errors to 1 sigma
fitresult_int = confint(fitresult, erf(1 / sqrt(2)));
fitresult_int = diff(fitresult_int,1)/2;

% save local parameters
a85 = fitresult.a85; %[V]
a85_err = fitresult_int(1);%[V]
a87 = fitresult.a87; %[V]
a87_err = fitresult_int(2);%[V]
G85 = fitresult.g85; %[Hz]
G87 = fitresult.g87; %[Hz]
G85_err = fitresult_int(6);    %[Hz]
G87_err = fitresult_int(7);    %[Hz]
T2_85 = 1 / G85; %[sec]
T2_85_err = G85_err / G85 .^ 2; % 
T2_87 = 1 / G87; %[sec]
T2_87_err = G87_err / G87 .^ 2; % 
f85 = fitresult.f85; %[Hz]
f85_err = fitresult_int(4);    %[Hz]
f87 = fitresult.f87; %[Hz]
f87_err = fitresult_int(5);    %[Hz]
phi85 = fitresult.phi85; %[rad / sec]
phi85_err = fitresult_int(8);  %[rad / sec]
phi87 = fitresult.phi87; %[rad / sec]
phi87_err = fitresult_int(9); %[rad / sec]



% Plot
figure(1);
subplot(211)
h = plot(fitresult, xData, yData);
legend('V(t)', 'fit', 'Location', 'Best' );
xlabel('t (sec)')
ylabel('v(t) (Volt)')
title(['T_2 Rb85 = ' num2str(T2_85*1e6, '%.4f') '+/-' num2str(T2_85_err*1e6, '%.4f') 'us, '...
    'T_2 Rb87 = ' num2str(T2_87*1e6, '%.4f') '+/-' num2str(T2_87_err*1e6, '%.4f') 'us'])
subplot(212)
findpeaks(P1, f, 'Annotate', 'extents', 'MinPeakProminence', 0.1 * max(P1(:)));
xlim([0 2 * max(locs)])
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')

% save data to meas structure
meas.FID_Bz =  f85 * 2 * pi / abs(c.gRb85); % Gauss
meas.FID_T2_85 = T2_85; % sec
meas.FID_T2_85_err = T2_85_err; % sec
meas.FID_a85 = a85;
meas.FID_a85_err = a85_err;
meas.FID_T2_87 = T2_87; % sec
meas.FID_T2_87_err = T2_87_err; % sec
meas.FID_a87 = a87;
meas.FID_a87_err = a87_err;
meas.FID_f85 = f85;
meas.FID_f85_err = f85_err;
meas.FID_f87 = f87;
meas.FID_f87_err = f87_err;
meas.FID_phi85 = phi85;
meas.FID_phi85_err = phi85_err;
meas.FID_phi87 = phi87;
meas.FID_phi87_err = phi87_err;

    

logger.log('Alkali FID measurement done.')
