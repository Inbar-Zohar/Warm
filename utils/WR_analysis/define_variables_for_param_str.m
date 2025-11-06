function [names, units] = define_variables_for_param_str(type, prefactors)

if ~exist('prefactors','var') || isempty(prefactors)
    prefactors = [1e-3 1 1 1 1e3];
end

common_names = {'ARW','BI','RRW','RR'};
common_units = {'deg/\surdhr','deg/hr','deg/hr^{3/2}','deg/hr^2'};

switch type
    case 'allan'
        names = [{'AWN'},common_names(:)'];
        units = [{'deg'},common_units(:)'];
    case 'srPSD'
        names = [{'AWNd'},common_names(:)'];
        units = [{'deg/\surdHz'},common_units(:)'];
            case 'both'
        names = [{'AWN','AWNd','NBW'},common_names(1:2),'BI_{eff}',common_names(3:4)];
        units = [{'deg','deg/\surdHz','Hz'},common_units(1:2),'deg/hr',common_units(3:4)];
    otherwise
        error('Unknown type %s',type);
end
