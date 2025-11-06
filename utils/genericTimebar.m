function genericTimebar(num_done, num_total, frac_done, loopTic, msg, no_backspace)
	% num_done and num_total can be arrays for nested loops; order them
	% from outermost loop to innermost [outer_done, inner_done] etc.

	persistent last_msg_length;
	persistent last_msg_time;
	if isempty(last_msg_length)
		last_msg_length = 0;
	end
	if isempty(last_msg_time)
		last_msg_time = 0;
	end
	
	if nargin == 0
		return;
	end
	
	if nargin < 6
		no_backspace = 0;
	end
	
	if length(num_done) > 1
		real_num_done = num_done(1);
		for ii = 2:length(num_done)
			real_num_done = (real_num_done-1)*num_total(ii) + num_done(ii);
		end
		num_done = real_num_done;
		num_total = prod(num_total);
	end
	
	if num_done == 1 || no_backspace == 1
		last_msg_length = 0;
	end
	
	if isempty(frac_done)
		frac_done = num_done / num_total;
	end
	
	if nargin < 5
		msg = 'Done run';
	end
	
	thisToc = toc(loopTic);
	timenow = clock;  
	timenow = timenow(4)*60*60 + timenow(5)*60 + timenow(6);
	
	if timenow > last_msg_time && (timenow - last_msg_time) < 0.5 && frac_done ~= 1
		return;
	end
	
	time_left = (thisToc/frac_done)*(1-frac_done);
	eta_time = timenow + time_left;
	if isequal(num_total, [])
		disp_msg = sprintf('%s %.0f percent (elapsed %s, remaining %s, eta: %s)...\n', msg, 100*frac_done, timeToStr(thisToc), timeToStr(time_left), timeToStr(eta_time));
	else
		disp_msg = sprintf('%s %d/%d (elapsed %s, remaining %s, eta: %s)...\n', msg, num_done, num_total, timeToStr(thisToc), timeToStr(time_left), timeToStr(eta_time));
	end
	fprintf([repmat(sprintf('\b'), 1, last_msg_length), disp_msg]);
	last_msg_length = length(disp_msg);
	last_msg_time = timenow;
	
	if frac_done == 1
		last_msg_length = 0;
	end
end