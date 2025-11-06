function out_str = create_variable_str(name, val, err, prefactor, unit)

if ~isempty(name)
    name_str = [name, ' = '];
else
    name_str = '';
end

if isfinite(err) 
    val_str = numerr2str(val/prefactor, err/prefactor);
else
    val_str = num2str(val/prefactor,2);
end

switch prefactor
    case 1e-3
        unit_str = [' m' unit];
    case 1e3
        unit_str = [' k' unit];
    otherwise
        unit_str = [' ' unit];
end

out_str = [name_str val_str unit_str];
