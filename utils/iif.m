function res = iif(cond, res_true, res_false)
	if cond
		res = res_true;
	else
		res = res_false;
	end
end