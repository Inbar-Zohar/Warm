function struct2var(s)
    fields = fieldnames(s);
    for i = 1:numel(fields)
        assignin('caller', fields{i}, s.(fields{i}));
    end
end
