
% uncomment to plot waveform 
% meas=create_fig(meas,2); plot(meas.FID_Rb_t, meas.FID_Rb_v); xlabel('time'); ylabel('amplitude [volts]'); hold on

% find measured signal start index
% [~, ind1] = max(abs(meas.FID_Rb_v));
ind1 = find(meas.FID_Rb_t > 2.5e-6,1);

% % plot truncated signal - uncomment for plotting
% meas=create_fig(meas); plot(meas.FID_Rb_t(ind1:end), meas.FID_Rb_v(ind1:end)); hold on

% set variables for fitting process
t4a = meas.FID_Rb_t(ind1:end);
v4a = meas.FID_Rb_v(ind1:end);

[f, P1] = fft1(t4a, v4a);


try
    % Find peaks in FFT
    [pks, locs, w, p] =  findpeaks(P1, f, 'Annotate', 'extents',...
        'MinPeakProminence', 0.1 * max(P1(:)), 'SortStr', 'descend');

    inds2save = locs < 300e3;
    pks = pks(inds2save);
    locs = locs(inds2save);
    w = w(inds2save);
    p = p(inds2save);



    %If only the 85Rb peak was detected, estimate the 87Rb peak frequency
    if numel(locs) < 2
        gF_Rb87 = 0.7; %MHz/G
        gF_Rb85 = 0.47; %MHz/G    
        locs(2) = locs(1)*gF_Rb87/gF_Rb85;
    elseif numel(locs) > 2

    end

    % get envelope - T2
    %finding np2 
    %"The envelopes are determined using spline interpolation over local maxima separated by at least np samples".
    %the meaning is how many points there are between the peaks
    %this is how the envelope function can know how to "sense" the plot
    est_period = 1/mean([locs(1), locs(2)]);
    np2 = round((est_period/diff(t4a(1:2)))*2);
    [up2, lo2] = envelope(v4a, np2, 'peak');

    % estimate T2 by finding the time where the envelope decreases to half
    ind2 = find(up2 - up2(end) < (up2(1) - up2(end)) / 4, 1);
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
    meas=create_fig(meas,1);
    subplot(211)
    h = plot(fitresult, xData, yData);
    legend('V(t)', 'fit', 'Location', 'Best' );
    xlabel('t (sec)')
    ylabel('v(t) (Volt)')
    % title({['T_2 Rb85 = ' num2str(T2_85*1e6, '%.4f') '\pm' num2str(T2_85_err*1e6, '%.4f') ' \mus, '...
    %     'T_2 Rb87 = ' num2str(T2_87*1e6, '%.4f') '\pm' num2str(T2_87_err*1e6, '%.4f') ' \mus'],...
    %     ['a_{85} = ' num2str(a85, '%.4f') '\pm' num2str(a85_err, '%.4f') ' V, ' ...
    %      'a_{87} = ' num2str(a87, '%.4f') '\pm' num2str(a87_err, '%.4f') ' V']})
    title({char(meas.DateTime_var),...
        ['T_2 Rb85 = ' numerr2str(T2_85*1e6,T2_85_err*1e6) ' \mus, '...
        'T_2 Rb87 = ' numerr2str(T2_87*1e6,T2_87_err*1e6) ' \mus'],...
        ['a_{85} = ' numerr2str(a85,a85_err) ' V, ' ...
         'a_{87} = ' numerr2str(a87,a87_err) ' V']})

    grid on
    subplot(212)
    findpeaks(P1, f, 'Annotate', 'extents', 'MinPeakProminence', 0.1 * max(P1(:)));
    xlim([0 2 * max(locs)])
    % title(['f_{85} = ' num2str(f85*1e-3, '%.4f') '\pm' num2str(f85_err*1e-3, '%.4f') 'kHz, '...
    %     'f_{87} = ' num2str(f87*1e-3, '%.4f') '\pm' num2str(f87_err*1e-3, '%.4f') 'kHz'])
    title(['f_{85} = ' numerr2str(f85*1e-3,f85_err*1e-3), ' kHz, '...
        'f_{87} = ' numerr2str(f87*1e-3,f87_err*1e-3), ' kHz'])

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
    meas.FID_FOM85 = a85^2*T2_85;
    meas.FID_FOM85_err = meas.FID_FOM85 .* sqrt((2*a85_err/a85)^2 + (T2_85_err/T2_85)^2);
    meas.FID_FOM87 = a87^2*T2_87;
    meas.FID_FOM87_err = meas.FID_FOM87 .* sqrt((2*a87_err/a87)^2 + (T2_87_err/T2_87)^2);

catch e
    meas=create_fig(meas,1);
    subplot(211)
    h = plot(t4a, v4a,'.-');
    xlabel('t (sec)')
    ylabel('v(t) (Volt)')
    title(char(meas.DateTime_var))
    
    grid on
    subplot(212)
    findpeaks(P1, f, 'Annotate', 'extents', 'MinPeakProminence', 0.1 * max(P1(:)));
    xlim([0 4e5])
    xlabel('f (Hz)')
    ylabel('|P1(f)|')
    
    logger.log(['Failed fitting Alkali FID results.']);
    logger.log(['A problem occured in ' e.stack(1).name ' (line '  num2str(e.stack(1).line) ')']);
    logger.log(['The message was: ' e.message])
end
    


% logger.log('Alkali FID measurement done.')
