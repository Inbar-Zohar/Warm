function plot_timeseries_fft_psd(t_vec,val_vec,new_plot,do_log_x,movmean_param)

if nargin < 3 || isempty(new_plot)
    new_plot = true;
end
if nargin < 4 || isempty(do_log_x)
    do_log_x = 1;
end
if nargin < 5 || isempty(movmean_param)
    movmean_param = 1;
end


if new_plot
    if isnumeric(new_plot)
        myfig(new_plot);
    else
        myfig;
    end
    sp1 = mysubplot(3,1,1);
    xlabel('Time (sec)');
    ylabel('Values');
    sp2 = mysubplot(3,1,2);
    xlabel('Frequency (Hz)');
    ylabel('FFT')
    sp3 = mysubplot(3,1,3);
    xlabel('Frequency (Hz)');
    ylabel('srPSD')

else
    sp1 = subplot(3,1,1);
    sp2 = subplot(3,1,2);
    sp3 = subplot(3,1,3); 
end

fs = 1 / mean(diff(t_vec));
[f,P1,pks,locks,Y] = fft1(t_vec,val_vec,1);
[psd,fff]=pwelch(val_vec(2:end), [], [], [], fs);
srPSD = sqrt(psd);


%time domain
filt_val_vec = movmean(val_vec,movmean_param);
% plot(sp1,t_vec(1:movmean_param:end),val_vec(1:movmean_param:end));
plot(sp1,t_vec,filt_val_vec);


%fft
plot(sp2,f,P1)
if do_log_x
    sp2.XScale = 'log';
end
sp2.YScale = 'log';

%psd
plot(sp3,fff,srPSD)

if do_log_x
    sp3.XScale = 'log';
end
sp3.YScale = 'log';
