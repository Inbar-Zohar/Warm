function WPstruct = WP_to_text(prm)

%Display all the fields of the flattened struct
initstruct = struct();
WPstruct = flatten_struct(prm.WP, initstruct);
end


function outstruct = flatten_struct(instruct, initstruct, prefix)
    if nargin < 3, prefix = ''; end
    outstruct = initstruct;
    fns = fieldnames(instruct);
    for I=1:numel(fns)
        fval = instruct.(fns{I});
        if isstruct(fval)
            outstruct = flatten_struct(fval, outstruct, [prefix fns{I} '__']);
        else
            outstruct.([prefix fns{I}]) = fval;
        end
    end
end