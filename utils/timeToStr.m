function s = timeToStr(t,showSeconds)
	if nargin < 2
		showSeconds = 1;
	end
	
	if t < 60
		s = sprintf('%.0fs', t);
	elseif t < 3600
		s = sprintf('%.0fm:%02.0fs', floor(t/60), rem(t,60));
	elseif t < 86400
		if showSeconds
			s = sprintf('%.0fh:%02.0fm:%02.0fs', floor(t/3600), floor(rem(t,3600)/60), rem(rem(t,3600),60));
		else
			s = sprintf('%.0fh:%02.0fm', floor(t/3600), floor(rem(t,3600)/60));
		end
	else
		if showSeconds
			s = sprintf('%.0fd:%02.0fh:%02.0fm:%02.0fs', floor(t/86400), floor(rem(t,86400)/3600), floor(rem(t,3600)/60), rem(rem(t,3600),60));
		else
			s = sprintf('%.0fd:%02.0fh:%02.0fm', floor(t/86400), floor(rem(t,86400)/3600), floor(rem(t,3600)/60));
		end
	end
end