function [num_err_str] = num_and_err(num, err)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function gets a measured number and the measurement error and
% returns a string of the number as an int if the log10(error)>0, else it 
% returns a string of the number until the error digit. Also in returns 
% the str(ceil(error))
%
%   Examples:
%       126e1[4e1] = num_and_err(1258.467, 41.2)
%       1258.47[1e-2] = num_and_err(1258.467, 0.012)
%       1258[3e0] = num_and_err(1258.4, 2.923)
%       6e-5[3e+0] = num_and_err(0.00005786, 2.923)
%       5.8e-05[1e-6] = num_and_err(0.00005786, 0.000000923)
%
%   NOT WORKING PROPERLY... NEEDS FIXING
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isnan(err)
    errstr = 'nan';
    numstr = strcat(num2str(ceil(num / (10^floor(log10(num))))), 'e', num2str(ceil(log10(num))));
else
    numstr = num2str(num, sign(log10(num)) * ceil(abs(log10(num)) - sign(log10(err)) * ceil(abs(log10(err)))));
    err_sig_dig = round(err * 10^abs(round(log10(err))));
    err_dec_sign = '+';
    if sign(log10(err)) == -1
        err_dec_sign = '-';
    end
    errstr = strcat(num2str(err_sig_dig), 'e', err_dec_sign, num2str(round(abs(log10(err)))));
end

num_err_str = strcat(numstr, '[', errstr, ']');
end

