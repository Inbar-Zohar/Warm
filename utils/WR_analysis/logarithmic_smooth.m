function [x,y] = logarithmic_smooth(x_full, y_full, method, opts)

if startsWith(method,'ieee')
    l = length(y_full);
    max_full_power_of_2 = ceil(log2(l));
    if max_full_power_of_2 <= 5
        x = x_full;
        y = y_full;
    else
        this_inds = 1:32;
        x(this_inds) = x_full(this_inds);
        y(this_inds) = y_full(this_inds);
        df = diff(x_full(1:2))/2;
        for ind=6:max_full_power_of_2
            this_av_factor = 2^(ind-5);
            this_av_window = [this_av_factor/2-1 this_av_factor/2];
            if ind < max_full_power_of_2
                this_inds_orig = (2^(ind-1)+1):2^ind;
                this_inds_select = this_av_factor/2:this_av_factor:2^(ind-1);
            else
                this_inds_orig = (2^(ind-1)+1):l;
                this_inds_select = this_av_factor/2:this_av_factor:(l-2^(ind-1)-1);
            end
            this_inds = this_inds(end) + (1:length(this_inds_select));
            
            this_f = x_full(this_inds_orig);
            if strcmp(method,'ieee')
                this_srPSD_smoothed = movmean(y_full(this_inds_orig),this_av_window);
            elseif strcmp(method,'ieee-log')
                this_srPSD_smoothed = 10.^movmean(log10(y_full(this_inds_orig)),this_av_window);
            end
            x(this_inds) = this_f(this_inds_select)+df;
            y(this_inds) = this_srPSD_smoothed(this_inds_select);
        end
    end
elseif strcmp(method, 'interp_and_movmean')
    xhat_full = logspace(log10(x_full(1)),log10(x_full(end)), ...
        opts.smooth_movmean * opts.output_length + 1);
    filtered_y = exp(movmean(...
        interp1(log(x_full), log(y_full), log(xhat_full), opts.interp_method),...
        opts.smooth_movmean));
    x = xhat_full(1:opts.smooth_movmean:end);
    y = filtered_y(1:opts.smooth_movmean:end);   
else
    error('Unknown smooth method')
end