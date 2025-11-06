function ax = get_ax_handles(figNumber)
    % Initialize the output as empty
    ax = [];
    
    % Check if the figure exists
    if ~isvalid(figure(figNumber))
        error('Figure %d does not exist.', figNumber);
    end

    % Get all the axes handles in the figure
    fig = figure(figNumber);
    axesHandles = findall(fig, 'Type', 'axes');

    % Check if there are at least two axes
    if length(axesHandles) < 2
        return;
    end

    % Get positions of all axes
    positions = get(axesHandles, 'Position');
    
    % We only need the first two axes for comparison
    pos1 = positions{1};
    pos2 = positions{2};

    % Determine the order based on vertical position (y-coordinate)
    if pos1(2) > pos2(2) % Compare y-coordinates (second element in position)
        ax = axesHandles(1);
    elseif pos1(2) < pos2(2)
        ax = axesHandles([2 1]);
    else
        % If y-coordinates are the same, compare x-coordinates (first element in position)
        if pos1(1) < pos2(1)
            ax = axesHandles(1);
        else
            ax = axesHandles([2 1]);
        end
    end
end
