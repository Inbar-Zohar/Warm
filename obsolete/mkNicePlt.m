function mkNicePlt()
% The function runs over the contents of the figure and changes things
% to make it nice and clear

% Costa Precision
tmp = gcf;
for i=1:length(tmp.Children)
    ax = tmp.Children(i);
    if isa(ax,'matlab.graphics.axis.Axes') ||...
            isa(ax,'matlab.graphics.illustration.Legend')
   
        set(ax,'FontSize',18);
    end
end

end