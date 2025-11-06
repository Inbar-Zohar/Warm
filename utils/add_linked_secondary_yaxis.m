function add_linked_secondary_yaxis(varargin)
    % add_linked_secondary_yaxis Adds a secondary y-axis to an existing primary  
    % axis with predefined scaling, linking the axes when zooming manually 
    % or through the ylim command (like in matplotlib). 
    % For it to work properly, always leave the axes with the left
    % (primary) y axis active. 
    % 
    % add_linked_secondary_yaxis(transformFunc) creates a secondary y-axis to 
    %   the current axis handle, where transformFunc defines the 
    %   relationship between the primary and secondary y-axes (e.g., @(y) 2*y).
    % add_linked_secondary_yaxis(primaryAx, transformFunc) applies the secondary 
    %   y-axis to the axis handle defined in primaryAx.
    % add_linked_secondary_yaxis(...,color1,color2) applies color1 and color 2 
    %   to the primary and secondary axes. Defaults to both being black.    
    %
    % Usage:
    %   figure;
    %   x = linspace(0, 10, 100);
    %   y1 = sin(x);
    %   y2 = cos(x);
    %   plot(x, y1);  % Plot primary data
    %   add_linked_secondary_yaxis(@(y) 2*y);  % Add secondary axis with scaling factor 2
    %   yyaxis right  % Switch to secondary axis 
    %   plot(x, y2);  % Plot secondary data
    %   yyaxis left   % Not doing that breaks things

    not_default_ax = ~isa(varargin{1},'function_handle');
    if not_default_ax
        primaryAx = varargin{1};
    else   
        primaryAx = gca;    
    end
    transformFunc = varargin{1 + not_default_ax};
    if nargin >= (3 + not_default_ax)
        primaryColor = color_vec(varargin{2 + not_default_ax});
        secondaryColor = color_vec(varargin{3 + not_default_ax});
    else
        primaryColor = color_vec('k');
        secondaryColor = color_vec('k');
    end
        
    
    % Switch to the primary y-axis if necessary
    yyaxis(primaryAx, 'left');

    % Create secondary y-axis
    yyaxis(primaryAx, 'right');
    yyaxis(primaryAx, 'left');

    primaryAx.YAxis(1).Color = primaryColor;
    primaryAx.YAxis(2).Color = secondaryColor;
    primaryAx.YAxis(2).Scale = primaryAx.YAxis(1).Scale;
    
    % Set up listener for zoom/pan updates on the primary y-axis
    addlistener(primaryAx, 'YLim', 'PostSet', @(src, event) updateSecondaryAxis(primaryAx, transformFunc));
    addlistener(primaryAx, 'XLim', 'PostSet', @(src, event) updateSecondaryAxis(primaryAx, transformFunc));
    zoomObj = zoom(ancestor(primaryAx,'figure'));
    set(zoomObj, 'ActionPostCallback', @(~, ~) updateSecondaryAxis(primaryAx, transformFunc));

    % Ensure the secondary axis is updated immediately
    updateSecondaryAxis(primaryAx, transformFunc);

end

function updateSecondaryAxis(primaryAx, transformFunc)
    % UPDATESECONDARYAXIS Updates the secondary y-axis based on the primary
    % axis limits and a transformation function
    
    % Switch to the primary axis
    yyaxis(primaryAx, 'left');
    yLimitsPrimary = ylim(primaryAx);  % Get primary y-axis limits

    % Apply the transformation to get secondary y-axis limits
    yLimitsSecondary = transformFunc(yLimitsPrimary);

    % Switch to the secondary axis and set the limits
    yyaxis(primaryAx, 'right');
    ylim(yLimitsSecondary);  % Update secondary y-axis limits
    yyaxis(primaryAx, 'left');

end
