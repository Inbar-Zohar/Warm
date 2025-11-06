function set_ax_lims(type, source, x_raw, y_raw, y_err_raw, ax, present_in_deg_hr_sqrtHz)

if ~exist('ax','var') || isempty(ax); ax=gca; else; axes(ax); end

if ~exist('y_err_raw','var') || isempty(y_err_raw); y_err_raw = zeros(size(y_raw)); end
if ~exist('present_in_deg_hr_sqrtHz','var') || isempty(present_in_deg_hr_sqrtHz)
        present_in_deg_hr_sqrtHz=1; 
end

if present_in_deg_hr_sqrtHz && strcmp(type,'srPSD')
    factor = 60;
else
    factor = 1;
end

switch source
    case 'data'
        x = x_raw;
        y = factor*y_raw;
        y_err = factor*y_err_raw;
    case 'model'
        switch type
            case 'allan'
                func = @noise_model_allan_func;
            case 'srPSD'
                func = @noise_model_srPSD_func;
            otherwise
                error('Unknown type %s',type);
        end
        x_range = x_raw;
        noise_term_vals = y_raw;
        
        which_terms = isfinite(noise_term_vals);
        x = logspace(log10(x_range(1)),log10(x_range(2)),100);
        y = factor*func(x,noise_term_vals(which_terms),which_terms);
        
        y_err = zeros(size(y));
    otherwise
        error('unknwon source %s', source)
end

switch type
    case 'allan'
        [miny, minyind] = min(y(1:floor(end*0.9)));
        minyval = miny - y_err(minyind);
    case 'srPSD'
        minyval = min(y);
end
maxyval = max(y);
maxxval = max(x);
minxval = min(x);

xlim(extend_range(minxval,maxxval,1));

yl = extend_range(minyval,maxyval,2);
ylim(yl)
% if log10(yl(2) / yl(1)) > 1 &&  log10(yl(2) / yl(1)) <= 4
%     set(ax,'YTick', 10.^[floor(log10(yl(1))) : 1 : ceil(log10(yl(2)))]);
% elseif log10(yl(2) / yl(1)) > 5
%     set(ax,'YTick', 10.^[2*floor(log10(yl(1))/2) : 2 : 2*ceil(log10(yl(2))/2)]);
% end

end

