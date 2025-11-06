function substitute_if_does_not_exist(variable_name,value)
%This function sets the variable "variabl_name" to "value" only if it does
%not exist yet
if evalin('base',['~exist(''' variable_name ''')'])
    if ischar(value)
        evalin('base',[variable_name '=''' value ''';'])
    else
        evalin('base',[variable_name '=' num2str(value) ';'])
    end
elseif evalin('base',['isnan(' variable_name ')'])
    if ischar(value)
        evalin('base',[variable_name '=''' value ''';'])
    else
        evalin('base',[variable_name '=' num2str(value) ';'])
    end
end
end

