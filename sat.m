function output = saturation(input,low,up)

if input > up
    output = up;
elseif input < low
    output = low;
else
    output = input;    
end