function meas = create_fig(meas,fig)

if nargin < 2
    fig_handle = figure;
else
    fig_handle = figure(fig);
end
fig_num = fig_handle.Number;

clf;
hold all;
box on;
grid on;

if ~isfield(meas,'figs')
    meas.figs = [fig_handle];
else
    meas.figs(end+1) = fig_handle;
end