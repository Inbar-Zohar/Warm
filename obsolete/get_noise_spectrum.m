function [] = get_noise_spectrum(data, fs, names)
% get_noise_spectrum computes the PSD of the data and plot it in a single
% figure in units of volt/srHz
if ~exist('names', 'var') || isempty(names) || isnumeric(names) || islogical(names)
    names = string(min(size(data)));
end

data_shape = size(data);

% make the data appear in the second dimension
if data_shape(2) < data_shape(1)
    data = data';
    data_shape = size(data);
end
n_signals = data_shape(1);

figure();
title('Power Spectral Desnsity of Data', 'fontweight', 'bold', 'FontSize', 40, 'FontName', 'Cambria Math');

for i = 1:n_signals

    % compute the PSD of data
    dv = data(i, :);
    [PSD_dv, f] = pwelch(dv - mean(dv), [], [], [], fs);

    % plot
    names(i) = regexprep(names(i), '_', ' ');
    loglog(f, sqrt(PSD_dv), 'displayname', names(i));
    if i == 1
        hold on
    end
end
xlabel('f(Hz)');
ylabel('$\sqrt{PSD} (\frac{V}{\sqrt{Hz}})$', 'interpreter', 'latex');
grid on
legend show

end

