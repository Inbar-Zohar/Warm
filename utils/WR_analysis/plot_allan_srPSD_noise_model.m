function plot_handle = plot_allan_srPSD_noise_model(x_range, noise_term_vals, varargin)
% Plots either Allan deviation or srPSD calculated curve. 
% Returns the handle of the model curve plot. 
%
% plot_allan_srPSD_noise_model(x_range, noise_term_vals, type) plots the 
%       calculated noise curve.
%
% plot_allan_srPSD_noise_model(x_range, noise_term_vals, noise_term_errs, type) 
%       plots also the calculated curve upper bound curve.
%
% Input parameters:
%   x_range: The range over which to plot the model curve. Either tau in seconds
%       if plotting Allan, or f in Hz if plotting srPSD.
%   noise_term_vals: A vector of the noise amplitudes, written as 
%       [AWN/d ARW BI RRW RR].
%   noise_term_errs: Uncertainty of noise amplitudes. 
%   type: a string argument equals either to 'allan' or 'srPSD'.
%
% Name/Value pairs:
%   ax: axis information to plot data on. Empty for current axes (gca), or 
%       an axes object. Defaults to empty.
%   plot_all: whether to plot the contribution of individual terms over the
%       range where they are relevant. Defaults to true.
%   allowed_contr_factor: The factor from the overall curve for which we 
%       consider a contribution "relevant" (if plot_all is true). 
%       Defaults to 0.25.
%   unit_prefactors: The prefactors to be used by the units of the noise
%       terms (1e-3 for m, 1e3 for k, etc.). Defaults to [1e-3 1 1 1 1e3].
%   N_points: Number of points for model curve plot. Defaults to 500. 
%   lw: linewidth for the model curve. Defaults to 2.  
%   present_in_deg_hr_sqrtHz: whether to plot srPSD data in deg/hr/srHz
%       instead of deg/sthr. Defaults to true. 
%   show_title: whether to show model parameters in figure title. 
%       Defaults to true. 
%   plot_colors: Color information for the model plot. Empty for figure's
%       default, a 1x3 rgb vector, or a 6x3 rbg matrix. The latter is used  
%       to define different colors for the different terms (if plot_all is
%       true). Defaults to empty. 


p = inputParser;
addRequired(p,'x_range',@(x) isnumeric(x) && (length(x) == 2));
addRequired(p,'noise_term_vals',@(x) isnumeric(x) && (length(x) == 5));

checkString = @(x) (isstring(x) || ischar(x)) && any(strcmp({'allan','srPSD'},x));
if isempty(varargin) || ...
        ((length(varargin) >= 1 && ~checkString(varargin{1})) && ...
        (length(varargin) >= 2 && ~checkString(varargin{2})))
    error('Type not defined')
elseif checkString(varargin{1}) 
    noise_term_errs = nan*noise_term_vals;
    type = varargin{1};
    varargin = varargin(2:end);
elseif checkString(varargin{2}) 
    noise_term_errs = varargin{1};
    type = varargin{2};
    varargin = varargin(3:end);
end


default_Npoints = 500;
default_unit_prefactors = [1e-3 1 1 1 1e3];

isAxesObject = @(x) isempty(x) || isa(x, 'matlab.graphics.axis.Axes');
addParameter(p,'ax',[],@(x)isAxesObject(x));

generalized_bool = @(x) (islogical(x) || isnumeric(x)) && isscalar(x);
addParameter(p,'plot_all',true,@(x) generalized_bool(x));
addParameter(p,'unit_prefactors',default_unit_prefactors,@(x) isnumeric(x) && (length(x) == 5));
addParameter(p,'allowed_contr_factor',0.25,@(x) isnumeric(x) && isscalar(x));
addParameter(p,'present_in_deg_hr_sqrtHz',true,@(x)generalized_bool(x));
addParameter(p,'show_title',true,@(x)generalized_bool(x));
addParameter(p,'N_points',default_Npoints,@(x) isnumeric(x) && (isempty(x) || isscalar(x)));
addParameter(p,'lw',2,@(x) isnumeric(x) && isscalar(x));

check_color_size = @(x) isempty(x) || ...
    (isnumeric(x) && (size(x,2) == 3) && ismember(size(x,1),[1 6]));
addParameter(p,'plot_colors',[],@(x) check_color_size(x));

parse(p, x_range, noise_term_vals, varargin{:});

struct2var(p.Results);

if ~exist('ax','var') || isempty(ax); ax=gca; else; axes(ax); end
if isempty(N_points); N_points=default_Npoints; end
if isempty(unit_prefactors); unit_prefactors=default_unit_prefactors; end

% if ~exist('noise_term_errs','var') || isempty(noise_term_errs)
%     noise_term_errs=nan(size(noise_term_vals)); 
% end
% 
% if ~exist('ax','var') || isempty(ax); ax=gca; else; axes(ax); end
% 
% if ~exist('plot_all','var') || isempty(plot_all); plot_all=1; end
% 
% if ~exist('allowed_contr_factor','var') || isempty(allowed_contr_factor); ...
%         allowed_contr_factor=0.25; 
% end
% if ~exist('lw','var') || isempty(lw); lw=1.5; end
% if ~exist('present_in_deg_hr_sqrtHz','var') || isempty(present_in_deg_hr_sqrtHz)
%         present_in_deg_hr_sqrtHz=1; 
% end
% if ~exist('show_title','var') || isempty(show_title); show_title=0; end

if present_in_deg_hr_sqrtHz && strcmp(type,'srPSD')
    factor = 60;
else
    factor = 1;
end

tf = ishold;
hold all;
if ~exist('plot_colors','var') || isempty(plot_colors)
    forced_color = 0;
    if tf
        color_ind = get(gca, 'ColorOrderIndex');
    else
        color_ind = 1;
    end
    plot_colors = linecolor(color_ind);
else
    forced_color = 1;
end

plot_color = plot_colors(1,:);
plot_individual_lines = (size(plot_colors,1) == 6) & size(noise_term_vals,1)==1;
if ~plot_individual_lines; plot_colors = plot_color; end
    
switch type
    case 'allan'
        func = @noise_model_allan_func;
    case 'srPSD'
        func = @noise_model_srPSD_func;
    otherwise
        error('Unknown type %s',type);
end


which_terms = isfinite(noise_term_vals);
x = logspace(log10(x_range(1)),log10(x_range(2)),N_points);
total_y = func(x,noise_term_vals(which_terms),which_terms);
if ~plot_individual_lines
    legend_str = create_legend_str(noise_term_vals, noise_term_errs, type, unit_prefactors);
else
    legend_str = 'Model';
end
plot_handle = plot(x,factor*total_y,'-','LineWidth',lw,'color',plot_color,'DisplayName',legend_str);

[names, units] = define_variables_for_param_str(type);

if plot_all
    for ind_term=find(which_terms)
        this_y_term = func(x,noise_term_vals(ind_term),(1:5) == ind_term);
        if plot_individual_lines
            this_color = plot_colors(ind_term+1,:);
            hand_vis = 'on';
            this_lw = lw;
            if show_title
                this_disp_name = names{ind_term};
            else
                this_disp_name = create_variable_str(names{ind_term}, ...
                    noise_term_vals(ind_term), noise_term_errs(ind_term), ...
                    unit_prefactors(ind_term), units{ind_term});
            end
        else
            this_color = plot_color;
            hand_vis = 'off';
            this_lw = lw-1;
            this_disp_name = '';
        end
        if any(this_y_term ./ total_y > allowed_contr_factor)
            this_inds = this_y_term ./ total_y > allowed_contr_factor;
            this_x = x(this_inds);
            this_y = factor*this_y_term(this_inds);
        else
            this_x = nan; this_y = nan;
        end
        plot(this_x,this_y,...
            '--','color',this_color,'LineWidth',this_lw,...
            'HandleVisibility',hand_vis,'DisplayName',this_disp_name)
    end
    
    if ~isempty(noise_term_errs)
        upper_y = func(...
            x,noise_term_vals(which_terms)+noise_term_errs(which_terms),which_terms);
        plot(x,factor*upper_y,':','LineWidth',lw, 'color',plot_color,'HandleVisibility','off');
    end
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

if show_title
    title_str = create_title_str(noise_term_vals, noise_term_errs, type, unit_prefactors);
    title(title_str);
end

if ~forced_color
    set(gca, 'ColorOrderIndex',color_ind+1);
end
if ~tf
    hold off;
end