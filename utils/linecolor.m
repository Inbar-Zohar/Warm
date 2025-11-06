function rgb = linecolor(ind,cmap)
	if nargin < 2
		cmap = lines(7);
		rgb = cmap(mod(ind-1,7)+1,:);
	else
		rgb = cmap(ind,:);
	end
end