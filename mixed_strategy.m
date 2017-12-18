function visit = mixed_strategy(p_aggre,na)


% agent "row" needs to choose an agent to visit.
% this is determined from the probability matrix.

      
% we pick a random number between 0 and 1
x = rand(1);
      
% here we have to make a decision
for col = 1 : na    % for each one of the other agents
          
     if col == 1   &   x <= p_aggre(1,col)
              
         % then he will visit col = 1
         visit = col;
              
     elseif col > 1   &   col <= na-1   &   x > p_aggre(1,col-1)  ...
             &  x <= p_aggre(1,col)
              
         visit = col;
         
     elseif col == na-1  &  x > p_aggre(1,col)
         
         visit = na;
              
     end
end

