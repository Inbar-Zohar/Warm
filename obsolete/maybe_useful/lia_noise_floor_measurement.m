logger.log('LIA noise floor measuremnet.')
disp('[+] LIA noise floor measuremnet.');

measurements = 7;
senetivities = zeros(1, measurements);
snrs = zeros(1, measurements);
white_noise_level = zeros(1, measurements);
input_ranges = zeros(1, measurements);

% parameters
channels = [1];
lia_time_constant = input('What is the LIA time constant ?');
enbw = 5 / 64 / lia_time_constant;
full_scale = 10;
under_sample = 100000;

figure();
for i = 1:measurements
    sensitivity = input('What is the LIA sensitivity in Vrms ?');
    input_range = input('What is the unput range ?');
    gain = full_scale / sensitivity;
    
    % set scope 
    Tscale = 2;
    Tdelay = 0;
    scope1.setTscale(Tscale);
    scope1.setTdelay(Tdelay)
    scope1.setTref('CENTer');
    scope1.TrigSource('EXT');
    scope1.Single;
    pause(5 * Tscale);

    % start measurement
    scope1.ForceTrig;
    scope1.readyToRead;


    % read signal from scope
    [t_, v_] = scope1.Read(channels);
    t = t_(1:round(length(t_) / under_sample):end);
    v = v_(:, 1:round(length(v_) / under_sample):end);
    v = v / gain;

    % release scope
    scope1.setTref('LEFT'); 
    scope1.TrigSource('LINE'); 
    scope1.Run;
    
    fs = 1 / mean(diff(t));
    [PSD_dv, f] = pwelch(v - mean(v), [], [], [], fs);

    % plot
    loglog(f, sqrt(PSD_dv), 'displayname', strcat('sensitivity: ', num2str(sensitivity), ' Vrms'));
    if i == 1
        hold on
        title('Power Spectral Desnsity of Data', 'fontweight', 'bold',...
        'FontSize', 20, 'FontName', 'Cambria Math');
    end
    white_noise_spectrum = f(f < enbw);
    white_noise_spectrum = white_noise_spectrum(0 < white_noise_spectrum);
    white_noise_mean = mean(sqrt(PSD_dv(1:length(white_noise_spectrum))));
    white_noise_std = std(sqrt(PSD_dv(1:length(white_noise_spectrum))));
    truncated_psd = PSD_dv(sqrt(PSD_dv) < white_noise_mean + 2 * white_noise_std);
    white_noise = mean(sqrt(truncated_psd(1:length(white_noise_spectrum))));
        
    loglog(white_noise_spectrum, white_noise * ones(length(white_noise_spectrum), 1),...
        'displayname', strcat('sensitivity: ', num2str(sensitivity),... 
    ' Vrms - ', 'wn: ', num2str(white_noise), ' V/srHz, '));

    snrs(i) = sensitivity / white_noise;
    senetivities(i) = sensitivity;
    input_ranges(i) = input_range;
end
xlabel('f(Hz)');
ylabel('$\sqrt{PSD} (\frac{V}{\sqrt{Hz}})$', 'interpreter', 'latex');
grid on
legend show


figure()
yyaxis left
loglog(senetivities, snrs, '-o')
hold on
loglog(senetivities, ones(1,length(senetivities)) * (1 / 400e-9))
hold off
title('LIA SNR vs LIA sensitivity', 'fontweight', 'bold',...
    'FontSize', 20, 'FontName', 'Cambria Math');
ylabel('SNR [srHz]')
yyaxis right
loglog(senetivities, senetivities ./ snrs, '-o')
xlabel('Sensitivity [Vrms]')
ylabel('Noise Floor [V/srHz]')
grid on
