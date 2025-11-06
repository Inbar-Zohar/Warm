function dt_str = datetime_str(dt)
% The function accepts datetime variable and returns date string 
    
    if nargin < 1
        dt = datetime;
    end
    formatOut = 'yyyy_mm_dd';
    dt_str = datestr(dt, formatOut);
end