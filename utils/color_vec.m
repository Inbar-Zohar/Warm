function colorVec = color_vec(inputColor)
    % GETCOLOR Converts color input to a 3-element RGB vector.
    %   Accepts both standard color strings ('k', 'r', etc.) and RGB vectors.
    
    % Define a map of standard color string abbreviations to RGB vectors
    colorMap = containers.Map({'k', 'r', 'g', 'b', 'c', 'm', 'y', 'w'}, ...
                              {[0 0 0], [1 0 0], [0 1 0], [0 0 1], [0 1 1], [1 0 1], [1 1 0], [1 1 1]});
    
    if ischar(inputColor)  % If input is a color string
        if colorMap.isKey(inputColor)
            colorVec = colorMap(inputColor);
        else
            error('Invalid color string. Use standard abbreviations like ''k'', ''r'', etc.');
        end
    elseif isnumeric(inputColor) && length(inputColor) == 3  % If input is an RGB vector
        colorVec = inputColor;
    else
        error('Input must be either a 3-element RGB vector or a standard color string.');
    end
end
