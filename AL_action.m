function [ P, action_out ] = AL_action( ASP, R, P, p_asp, c_asp, na, action )

% na            = No_actions(s,1)
% action        = action(s,1)
% P             = P{s,1}
% ASP           = ASP{s,1}
% R             = R{s,1}

 if R >= ASP | na == 1
        action_out = action;
 else
    %x=ASP{s,1}-R{s,1};
    % probability of selecting the previous action
    P(action,1) = sigmoid_Karandikar(ASP-R,p_asp,c_asp);
    % probability of selecting any other action
    for cc = 1 : na;
        if cc ~= action
            P(cc,1) = (1 - P(action,1))/(na-1);
        end
    end

    % Create aggregate probability vector
    P_aggre = zeros(na,1);
    sum_p = 0;  % sum of probability entries
    for j = 1 : na
        P_aggre(j,1) = P(j,1) + sum_p;
        sum_p = P_aggre(j,1);		% update of the sum
    end

    % Players first choose an action
    action_out = mixed_strategy(P_aggre',na);
 end
    
end

