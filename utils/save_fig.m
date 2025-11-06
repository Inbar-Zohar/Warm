function save_fig(figname, format, include_datetime)
%                 (string, string, boolean)
% The function saves the curent figure with figname plus data and time at a
% format 'png' or 'pdf'.  Defaults:
%                             figname = curr_fig   
%                             format = png
%                             include_datetime = 1  
%  Costa Precision™

if nargin<1
    figname = 'curr_fig';
end
if nargin<2
    format = 'png';
end
if nargin<3
    include_datetime = 1;
end

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
if include_datetime
    opt.FileName = [figname '_' datetime_str '.' format];
else
    opt.FileName = [figname '.' format];
end
print(gcf, ['-d' format], opt.FileName);
end