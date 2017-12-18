% Random selection of an element of a vector

function out = random_selection(x)

% x is assumed a column vector

[choices,xa] = size(x);

prob = zeros(choices,1); % define a probability vector

for i = 1 : choices
    
    prob(i) = (i-1)*1/choices;
    
end
prob(i+1)=1;

xaxa = rand(1);

for i = 1 : choices
    if xaxa >= prob(i) & xaxa < prob(i+1)
        selection = i;
    end
end

out = x(selection,1);