function out_str = numerr2str(v, dv)

max_val_for_two_digits = 1.94;

leading_digit = @(x) floor(log10(abs(x*1.05)));
leading_digit_val = @(x) x / (10^leading_digit(x));
round_leading_digit_val = @(x) round(leading_digit_val(x), ...
    single(leading_digit_val(x)< max_val_for_two_digits));

if dv == 0 || ~isfinite(dv)
    out_str = num2str(v);
else
    try
        dv = abs(dv);
        n_dv = leading_digit(dv);
        if leading_digit_val(dv) < max_val_for_two_digits
            n_dv = n_dv - 1;
        elseif leading_digit_val(dv) >= 9.5            
            n_dv = n_dv + 1;
        end
        if dv / (10^n_dv) < max_val_for_two_digits
            n_dv = n_dv - 1;
        end

        round_v = round(v,-n_dv);
        uncert = round_leading_digit_val(dv);

        if ( ((abs(round_v) >= 1e-2 && abs(round_v) < 1e3) || round_v == 0) ...
                && dv < 9.5  ) 
            if n_dv >= 0
                val_str = num2str(round_v);
            else
                eval(['val_str = sprintf("%.',num2str(-n_dv), 'f",v);']); 
            end

            if uncert < max_val_for_two_digits && n_dv ~= -1
                uncert = uncert*10;
            end

            suf_str = '';
        else
            n_exp = max(n_dv, leading_digit(round_v));
            if uncert >= 1 && uncert <= max_val_for_two_digits && n_dv >= leading_digit(round_v)
                n_exp = n_exp + 1;
            end

            eval(['val_str = sprintf("%.',num2str(n_exp-n_dv), 'f",v/10^n_exp);'])

            if uncert < max_val_for_two_digits && n_exp ~= n_dv + 1
                uncert = uncert*10;
            end            


            suf_str = sprintf(' x 10^{%i}',n_exp);
        end
        switch uncert
            case 1
                if val_str(end) == '0'
                    val_str = val_str(1:end-2);
                    uncert_str = '1';
                else
                    uncert_str = '1.0';
                end
            case 10 
                if val_str(end) == '0'
                    val_str = val_str(1:end-1);
                    uncert_str = '1';
                else
                    uncert_str = '10';
                end
            otherwise
                uncert_str = num2str(uncert);
        end
        out_str = sprintf('%s(%s)%s',val_str,uncert_str,suf_str);
    catch
        out_str = [num2str(v) '+/-' num2str(dv)]; 
        disp(['Something weird happened when trying to do numerr2str for ', out_str])
    end
        
end