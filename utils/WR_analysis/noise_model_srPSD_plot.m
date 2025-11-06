function srPSD_plot(f_data, srPSD_data, srPSD_err_data, frangefit, noise_term_vals, noise_term_errs, N_points, factor, ax)

if ~exist('ax','var') || isempty(ax)
    ax=gca;
else
    axes(ax);
end

if ~exist('factor','var') || isempty(factor)
    factor=0.25;
end

if ~exist('N_points','var') || isempty(N_points)
    N_points=500;
end

% if ~exist('noise_term_errs','var') || isempty(noise_term_errs)
    noise_term_errs = nan(size(noise_term_vals));
% end


tf = ishold;

plot(f_data, srPSD_data, '-', 'Color', linecolor(1),'LineWidth',2.5);
% patch('XData',[f_all;flip(f_all)], 'YData', [srPSD_all-srPSD_err_all;flip(srPSD_all+srPSD_err_all)],...
%         'EdgeColor', 'none', 'FaceColor', linecolor(1),'FaceAlpha',0.2);

which_terms = isfinite(noise_term_vals);
f = logspace(log10(frangefit(1)),log10(frangefit(2)),N_points);
total_srPSD = srPSD_func(f,noise_term_vals(which_terms),which_terms);
plot(f,total_srPSD,'-','LineWidth',1.5, 'color',linecolor(2));

for ind=find(which_terms)
    this_srPSD_term = srPSD_func(f,noise_term_vals(ind),(1:5) == ind);
    if any(this_srPSD_term ./ total_srPSD > factor)
        this_inds = this_srPSD_term ./ total_srPSD > factor;
        plot(f(this_inds),this_srPSD_term(this_inds),'--','color',linecolor(2),'LineWidth',0.5)
    end
end
% 
% upper_AD = srPSD_func(f,noise_term_vals(which_terms)+noise_term_errs(which_terms),which_terms);
% plot(Ts,upper_AD,':','LineWidth',1.5, 'color',linecolor(2));


title_str = {'',''};
names = {'AWNd','ARW','BI','RRW','RR'};
units = {'mdeg/\surdHz','deg/\surdhr','deg/hr','deg/hr^{3/2}','kdeg/hr^2'};
prefactors = [1e3 1 1 1 1e-3];
for ind = find(which_terms)
    prefactor = prefactors(ind);
    if isfinite(noise_term_errs(ind))
        this_str = sprintf('%s = %s %s', names{ind}, ...
            numerr2str(noise_term_vals(ind)*prefactor, noise_term_errs(ind)*prefactor), units{ind});
    else
        this_str = sprintf('%s = %s %s', names{ind}, ...
            num2str(noise_term_vals(ind)*prefactor,2), units{ind});
    end
    
    this_row = iif(ind <= 2, 1, 2);
    title_str{this_row} = [title_str{this_row} '; ', this_str];
end
for row=1:2
    if ~isempty(title_str{row})
        title_str{row} = title_str{row}(3:end);
    end
end

if ~isempty(srPSD_data)
%     [minsrPSD, minsrPSDind] = min(srPSD_all);
%     minval = minAD - AD_err_all(minADind);
    maxval = max(srPSD_data);
    minval = min(srPSD_data);
    f_for_lim = f_data;
else
    minval = min(total_srPSD);
    maxval = max(total_srPSD);
    f_for_lim = f;
end

pow10min = floor(log10(minval));
leaddigitmin = floor(minval/10.^pow10min);
if leaddigitmin < 3
    leaddigitmin = leaddigitmin*10;
    pow10min = pow10min-1;
end
ylmin = floor(leaddigitmin/3)*10.^pow10min;
pow10max = floor(log10(maxval));
leaddigitmax = floor(maxval/10.^pow10max);
ylmax = leaddigitmax*2*10.^pow10max;

xlmin = min(f_for_lim);
xlmax = max(f_for_lim);

ylim([ylmin ylmax])
xlim([xlmin xlmax])
title(title_str)
ylabel('\surdPSD_\Omega (deg/\surdhr)')
xlabel('f (Hz)')

if ~tf
    hold off;
end

