function [ action_out ] = BB_action ( action, action_base, s, epsilon_bb, No_actions, t )

if ( t > 1 )
    % we only select actions, after the first iteration, so that we
    % first define the baseline action and utility
    rn = rand(1);
    if (rn < epsilon_bb)
        % Random selection of an action
        % First create a uniform distribution
        % Create aggregate probability vector
        P_aggre = zeros(No_actions(s,1),1);
        sum_p = 0;  % sum of probability entries
        for j = 1 : No_actions(s,1)
            P_aggre(j,1) = 1/No_actions(s,1) + sum_p;
            sum_p = P_aggre(j,1);		% update of the sum
        end
        action_out = mixed_strategy(P_aggre',No_actions(s,1));
    else
        action_out = action_base(s,1);
    end
else
    action_out = action(s,1);
end

end

