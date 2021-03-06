function [ action_out ] = MB_action ( No_actions, action, action_base, s, lambda_mb, c_mb, t )

 if (t > 1)
    % we only select actions, after the first iteration, so that we
    % first define the baseline action and utility
    if ( mode(s,1) == 1 )
        % if this player is 'content'
        % we need to formulate a vector for selecting actions
        P_aggre = zeros(No_actions(s,1),1);
        sum_p = 0;  % sum of probability entries
        for j = 1 : No_actions(s,1)
            if ( j ~= action_base(s,1) )
                P_aggre(j,1) = lambda_mb^(c_mb)/(No_actions(s,1)-1) + sum_p;
                sum_p = P_aggre(j,1);		% update of the sum
            else
                P_aggre(j,1) = 1 - lambda_mb^(c_mb) + sum_p;
                sum_p = P_aggre(j,1);
            end
        end
    else
        % we need to formulate a vector for selecting actions
        P_aggre = zeros(No_actions(s,1),1);
        sum_p = 0;  % sum of probability entries
        for j = 1 : No_actions(s,1)
            P_aggre(j,1) = 1/No_actions(s,1) + sum_p;
            sum_p = P_aggre(j,1);		% update of the sum
        end
    end
    action_out = mixed_strategy(P_aggre',No_actions(s,1));
 else
     action_out = action(s,1);
 end

end

