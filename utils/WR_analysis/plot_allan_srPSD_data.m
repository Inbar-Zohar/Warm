function [plot_handle, patch_handle] = plot_allan_srPSD_data(x, y, varargin)
% Plots either Allan deviation or srPSD curve as a solid line, formatting
% the axes according to the specified data type. If uncertainty is supplied, 
% it is shown as a shaded area. 
% Returns the handle of the data plot and the uncertainty patch (latter is 
% empty if uncertainy not defined). 
%
% plot_allan_srPSD_data(x, y, type) plot vector y vs. vector x and formats
% figure according to type. 
%
% plot_allan_srPSD_data(x, y, y_err, type) adds symmetric uncertainty y_err
% around y.
%
% plot_allan_srPSD_data(x, y, y_err_neg, y_err_pos, type) uses y_err_neg 
% and y_err_pos for negative and positive uncertainty around y.
%
% Input parameters:
%   x, y: Can be either tau, AD in seconds, deg/hr for if plotting Allan, 
%       or f, srPSD in Hz, deg/srhr if plotting srPSD. 
%   y_err(/_neg/_pos): Uncertainty in y, same same and same units as y. 
%   type: a string argument equals either to 'allan' or 'srPSD'.
% 
% Name/Value pairs:
%   ax: axis information to plot data on. Empty for current axes (gca), or 
%       an axes object. Defaults to empty. 
%   plot_color: Color information for the data plot. Empty for figure's
%       default, or a 1x3 rgb vector. Defaults to empty. 
%   lw: linewidth for the data curve. Defaults to 2.5.  
%   patch_alpha: Face filling alpha value for uncertainty shaded area. 
%       Defaults to 0.2. 
%   legend_entry: A string describing the data, to be used in title/legend. 
%       Defaults to 'Data'.
%   present_in_deg_hr_sqrtHz: whether to plot srPSD data in deg/hr/srHz
%       instead of deg/sthr. Defaults to true. 

p = inputParser;
addRequired(p,'x',@isnumeric);
addRequired(p,'y',@isnumeric);

checkString = @(x) (isstring(x) || ischar(x)) && any(strcmp({'allan','srPSD'},x));
if isempty(varargin) || ...
        ((length(varargin) >= 1 && ~checkString(varargin{1})) && ...
        (length(varargin) >= 2 && ~checkString(varargin{2})) && ...
        (length(varargin) >= 3 && ~checkString(varargin{3})))
    error('Type not defined')
elseif checkString(varargin{1}) 
    y_err_neg = [];
    y_err_pos = [];
    type = varargin{1};
    varargin = varargin(2:end);
elseif checkString(varargin{2}) 
    y_err_neg = varargin{1};
    y_err_pos = varargin{1};
    type = varargin{2};
    varargin = varargin(3:end);
elseif checkString(varargin{3}) 
    y_err_neg = varargin{1};
    y_err_pos = varargin{2};
    type = varargin{3};
    varargin = varargin(4:end);    
end

isAxesObject = @(x) isempty(x) || isa(x, 'matlab.graphics.axis.Axes');
addParameter(p,'ax',[],@(x)isAxesObject(x));
addParameter(p,'plot_color',[],@(x) isempty(x) || (isnumeric(x) && (length(x) == 3)));
addParameter(p,'lw',2.5,@(x) isnumeric(x) && isscalar(x));
addParameter(p,'patch_alpha',0.2,@(x) isnumeric(x) && isscalar(x));
addParameter(p,'legend_entry','Data',@(x) isstring(x) || ischar(x));

generalized_bool = @(x) (islogical(x) || isnumeric(x)) && isscalar(x);
addParameter(p,'present_in_deg_hr_sqrtHz',true,@(x)generalized_bool(x));

parse(p, x, y, varargin{:});

struct2var(p.Results);

if ~exist('ax','var') || isempty(ax); ax=gca; else; axes(ax); end

% if isempty(y_err_pos); y_err_pos=y_err_neg; end

% if ~exist('lw','var') || isempty(lw); lw=2.5; end
% if ~exist('alpha','var') || isempty(alpha); alpha=0.2; end
% if ~exist('legend_entry','var'); legend_entry='Data'; end
% if ~exist('present_in_deg_hr_sqrtHz','var') || isempty(present_in_deg_hr_sqrtHz)
%         present_in_deg_hr_sqrtHz=1; 
% end


if present_in_deg_hr_sqrtHz && strcmp(type,'srPSD')
    factor = 60;
else
    factor = 1;
end

tf = ishold;
hold all;
if ~exist('plot_color','var') || isempty(plot_color)
    forced_color = 0;
    if tf
        color_ind = get(gca, 'ColorOrderIndex');
    else
        color_ind = 1;
    end
    plot_color = linecolor(color_ind);
else
    forced_color = 1;
end

if ~isempty(legend_entry)
    plot_handle = plot(x(:), factor*y(:), '-', 'Color',plot_color,'LineWidth',lw,'DisplayName',legend_entry);
else
    plot_handle = plot(x(:), factor*y(:), '-', 'Color',plot_color,'LineWidth',lw,'HandleVisibility','off');
end

patch_handle = [];
if ~isempty(y_err_neg)
    patch_handle = patch('XData',[x(:);flip(x(:))], 'YData', factor*[y(:)-y_err_neg(:);flip(y(:)+y_err_pos(:))],...
        'EdgeColor', 'none', 'FaceColor',plot_color,'FaceAlpha',patch_alpha,'HandleVisibility','off');
end

switch type
    case 'allan'
        xlabel('Integration time (s)')
        ylabel('\sigma_\Omega (deg/hr)')
    case 'srPSD'
        xlabel('f (Hz)')
        if present_in_deg_hr_sqrtHz
            ylabel('\surdPSD_\Omega (deg/hr/\surdHz)')            
        else
            ylabel('\surdPSD_\Omega (deg/\surdhr)')
        end
    otherwise
        error('Unknown type %s',type);
end
    

if ~forced_color
    set(gca, 'ColorOrderIndex',color_ind+1);
end
if ~tf
    hold off;
end