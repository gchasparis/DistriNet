function y=sigmoid_Karandikar(x,p,c)

% x = (1-p)/2*tanh(100*y+3) + (1+p)/2;

if x <= 0
    y=1;      
elseif x > 0
    y=max(p,-c*x+1); % here c=-2=-1/z_0, originally 0.2
end
