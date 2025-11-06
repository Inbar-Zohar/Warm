function sig = get_act_LIA_signal(raw_sig, ch_num, LIAconfig)

output_ch = LIAconfig.(['output_ch' num2str(ch_num)]);
if ismember(output_ch,{'X','Y','R'})
    if isfield(LIAconfig,['output_expand_factor_' lower(output_ch)])
        expand_factor = LIAconfig.(['output_expand_factor_' lower(output_ch)]);
    else
        expand_label = LIAconfig.(['output_expand_' lower(output_ch)]);
        switch expand_label
            case 'OFF'; expand_factor = 1;
            case 'x10'; expand_factor = 10;
            case 'x100'; expand_factor = 100;
        end
    end
    scaling = LIAconfig.sensitivity / 10 / expand_factor;

    if isfield(LIAconfig,['output_offset_act_' lower(output_ch)])
        offset_value = LIAconfig.(['output_offset_act_' lower(output_ch)]);
    else
        offset_stt = 1;
        if isfield(LIAconfig,['output_offset_' lower(output_ch) '_state'])
            offset_stt = LIAconfig.(['output_offset_' lower(output_ch) '_state']);
        end
        offset_percent = offset_stt * LIAconfig.(['output_offset_' lower(output_ch) '_value']);
        offset_value = LIAconfig.sensitivity * offset_percent / 100;
    end
else
    scaling = 180 / 10;
    offset_value = 0;
end

sig = raw_sig * scaling + offset_value;