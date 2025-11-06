function ax = plot_allan_srPSD_noise_wrapper(tau, AD, AD_err, f, srPSD, srPSD_err,...
    tau_range, f_range, noise_term_vals, noise_term_errs, varargin)

% A wrapper function for summary presentation of an analyzed experiment. It
% plots on two panels, one for Allan and one for srPSD, the experimental
% data and model (i.e., fit/s), with unceretainty if provided.
% Returns a vector of the two axes objects for the Allan and srPSD plots.
%
% Input parameters:
%   tau, AD, AD_err: Allan deviation data (integration time, ADEV, ADEV
%       uncertainty). tau in seconds, AD and AD_err in deg/hr.
%   f, srPSD, srPSD_err: srPSD data (frequency, srPSD, srPSD uncertainty).
%       f in Hz, srPSD in deg/srhr, srPSD_err a two-element vector
%       of relative bottom and top uncertainty.
%   tau_range, f_range: range over which the model should plotted
%   noise_term_vals, noise_term_errs: the noise model parameters and
%       uncertainties. Written as [AWN ARW BI RRW RR AWNd].
%
% Name/Value pairs:
%   fig: figure information to plot data on. Empty for a new figure, scalar
%       for creating a figure with this number, figure object, or a vector
%       of two axes objects. Defaults to empty.
%   reset_fig: whether to reset the figure if it exists. Defaults to true.
%   update_lims: whether to use optimize ax limits to data. Defaults to true.
%   plot_all: whether to plot the contribution of individual noise terms. Defaults to true.
%   show_title: whether to show model parameters in figure title. Defaults to true.
%   show_legend: whether to show a legend with model parameters as the
%       legend entry for the model. Defaults to false.
%   show_ngc: whether to plot the NGC NMRG performance curves. Defaults to true.
%   data_leg: A string describing the data, to be used in title/legend. Defaults to 'Data'.
%   unit_prefactors: The prefactors to be used by the units of the noise
%       terms (1e-3 for m, 1e3 for k, etc.). Defaults to [1e-3 1 1 1 1e3].
%   N_points: Number of points for model curve plot. Empty for using the default
%       value of the daughter function or an integer. Defaults to empty.
%   linecolor_data: Color information for the data plot. Empty for figure's
%       default, or a 1x3 rgb vector. Defaults to empty.
%   linecolor_model: Color information for the model plot. Empty for figure's
%       default, a 1x3 rgb vector for all model curves, or a 6x3 rgb matrix,
%       first two for the full model and the next for the 5 noise terms (assuming
%       plot_all is true). Defaults to empty.

p = inputParser;
addRequired(p,'tau',@isnumeric);
addRequired(p,'AD',@isnumeric);
addRequired(p,'AD_err',@isnumeric);
addRequired(p,'f',@isnumeric);
addRequired(p,'srPSD',@isnumeric);
addRequired(p,'srPSD_err',@(x) isempty(x) || (isnumeric(x) && (length(x) == 2)));
addRequired(p,'tau_range',@(x) isnumeric(x) && (length(x) == 2));
addRequired(p,'f_range',@(x) isnumeric(x) && (length(x) == 2));
% addRequired(p,'noise_term_vals',@(x) isnumeric(x) && (length(x) == 6));
addRequired(p,'noise_term_vals');
% addOptional(p,'noise_term_errs',nan*noise_term_vals,@(x) isnumeric(x) && (length(x) == 6));
addOptional(p,'noise_term_errs',nan*noise_term_vals);

validFig = @(fig) isempty(fig) || ...
    (isnumeric(fig) && isscalar(fig) && mod(fig,1)<1e-10 && fig > 0) || ...
    (length(fig) == 1 && isvalid(fig) && strcmp(get(fig, 'Type'), 'figure')) || ...
    (isvector(fig) && length(fig) >= 2 && all(isvalid(fig)) && all(strcmp(get(fig, 'Type'), 'axes')));
addParameter(p,'fig',[],@(x)validFig(x));

generalized_bool = @(x) (islogical(x) || isnumeric(x)) && isscalar(x);
addParameter(p,'reset_fig',true,@(x)generalized_bool(x));
addParameter(p,'update_lims',true,@(x)generalized_bool(x));
addParameter(p,'plot_all',true,@(x)generalized_bool(x));
addParameter(p,'show_title',true,@(x)generalized_bool(x));
addParameter(p,'show_legend',false,@(x)generalized_bool(x));
addParameter(p,'show_ngc',true,@(x)generalized_bool(x));

addParameter(p,'data_leg','Data',@(x) isstring(x) || ischar(x));
addParameter(p,'unit_prefactors',[1e-3 1e-3 1 1 1],@(x) isnumeric(x) && (length(x) == 5));

scalar_or_empty = @(x) (isnumeric(x) && isscalar(x)) || isempty(x);
addParameter(p,'N_points',[],@(x)scalar_or_empty(x));
addParameter(p,'linecolor_data',[],@isnumeric);
addParameter(p,'linecolor_model',[],@isnumeric);

parse(p, tau, AD, AD_err, f, srPSD, srPSD_err,...
    tau_range, f_range, noise_term_vals, noise_term_errs, varargin{:});

struct2var(p.Results);

if isempty(fig)
    fig_h=gcf;
    reset_fig = 1;
elseif isnumeric(fig) && isscalar(fig) || ...
        length(fig) == 1 && isvalid(fig) && strcmp(get(fig, 'Type'), 'figure')
    if isnumeric(fig) && isscalar(fig)
        fig_n = fig;
        figure(fig_n);
    else
        figure(fig);
        fig_n = fig.Number;
    end
    ax = get_ax_handles(fig_n);
    if isempty(ax); reset_fig = 1; end
elseif isvector(fig) && length(fig) >= 2 && ...
        all(isvalid(fig)) && all(strcmp(get(fig, 'Type'), 'axes'))
    ax = fig(1:2);
end

if reset_fig
    clf;
    try
        tiledlayout(2,1,'Padding','compact','TileSpacing','compact');
        for ind=1:2
            nexttile;
            ax(ind) = gca;
            hold all;
            grid on;
            box on;
            set(ax(ind),'XScale','log','YScale','log');
        end
    catch
        for ind=1:2
            ax(ind) = subplot(2,1,ind);
            hold all;
            grid on;
            box on;
            set(ax(ind),'XScale','log','YScale','log');
        end
    end
end

plot_individual_lines = (size(linecolor_model,1) == 6) & size(noise_term_vals,1)==1;
if plot_individual_lines; data_label = 'Data'; else; data_label = data_leg; end

axes(ax(1));
plot_allan_srPSD_data(tau, AD, AD_err, 'allan','plot_color',linecolor_data,'legend_entry',data_label);
for ind=1:size(noise_term_vals,1)
    plot_allan_srPSD_noise_model(tau_range, ...
        noise_term_vals(ind,1:5), noise_term_errs(ind,1:5),'allan', ...
        'plot_all', plot_all, 'unit_prefactors', unit_prefactors, 'N_points', N_points,...
        'plot_colors', linecolor_model);
end
if show_ngc
    ngc_noise_params = load_ngc_noise_params;
    plot_handle = plot_allan_srPSD_noise_model(tau_range, ngc_noise_params(1:5), 'allan', ...
        'plot_all', 0, 'show_title', 0, 'plot_colors', zeros(1,3));
    plot_handle.DisplayName = 'NGC';
    plot_handle.LineStyle = ':';
end
if show_title
    title(create_title_str(noise_term_vals(1,:), noise_term_errs(1,:), 'both', ...
        unit_prefactors, data_leg,tau_range))
end
if update_lims
    set_ax_lims('allan', 'data', tau, AD, AD_err);
end
if plot_individual_lines
    legend('location','southwest');
end

axes(ax(2));
plot_allan_srPSD_data(f, srPSD, srPSD.*srPSD_err(1), srPSD.*srPSD_err(2), 'srPSD',...
    'plot_color',linecolor_data,'legend_entry',data_leg);
for ind=1:size(noise_term_vals,1)
    plot_allan_srPSD_noise_model(f_range, ....
        noise_term_vals(ind,[6 2:5]), noise_term_errs(ind,[6 2:5]),'srPSD', ...
        'plot_all', plot_all, 'unit_prefactors', unit_prefactors, 'N_points', N_points,...
        'show_title', 0, 'plot_colors', linecolor_model);
end
if show_ngc
    plot_handle = plot_allan_srPSD_noise_model(f_range, ngc_noise_params([6 2:5]), 'srPSD', ...
        'plot_all', 0, 'show_title', 0, 'plot_colors', zeros(1,3));
    plot_handle.DisplayName = 'NGC';
    plot_handle.LineStyle = ':';
end
if update_lims
    set_ax_lims('srPSD', 'data', f, srPSD);
end
if show_legend
    legend('location','southoutside');
end
