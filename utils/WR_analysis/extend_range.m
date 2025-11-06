function out = extend_range(min_in,max_in,factor)

pow10min = floor(log10(min_in));
leaddigitmin = floor(min_in/10.^pow10min);
if leaddigitmin < factor
    leaddigitmin = leaddigitmin*10;
    pow10min = pow10min-1;
end
out(1) = floor(leaddigitmin/factor)*10.^pow10min;

pow10max = floor(log10(max_in));
leaddigitmax = ceil(max_in/10.^pow10max);
out(2) = leaddigitmax*factor*10.^pow10max;
end
