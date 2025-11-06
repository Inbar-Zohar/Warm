WPstruct = WP_to_text(prm);
fields = fieldnames(WPstruct);

fid = fopen('working_point.txt', 'wt');

for i = 1:numel(fields)
    fprintf(fid, '%s: ', fields{i}); % write field name
    value = WPstruct.(fields{i}); % get value of the field
    if ischar(value)
        fprintf(fid, '%s\n', value); % write string value
    elseif isnumeric(value)
        fprintf(fid, '%s\n', num2str(value)); % write numeric value
    else
        fprintf(fid, 'Unsupported data type\n');
    end
end

fclose(fid);
