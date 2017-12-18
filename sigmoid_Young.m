function y=sigmoid_Young(x,p)

% x = (1-p)/2*tanh(100*y+3) + (1+p)/2;
a = 0.1;
if ( a * x + 1/2 < 1 - p && a * x + 1/2 > p )
    y = a * x + 1/2;
elseif ( a * x + 1/2 >= 1-p )
    y = 1-p;
else
    y = p;
end
