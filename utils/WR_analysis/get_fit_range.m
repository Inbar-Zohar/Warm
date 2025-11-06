function [allan_min_ind,allan_max_ind,srPSD_min_ind,srPSD_max_ind] = ...
    get_fit_range(tau_s,f,est_bw)

allan_min_ind = find(tau_s>2/est_bw,1);
allan_max_ind = length(tau_s)-1;
srPSD_min_ind = 1;
srPSD_max_ind = find(f<est_bw/2.5,1,'last');
