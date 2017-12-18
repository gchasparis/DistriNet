function [action_out] = TE_action ( No_actions, action, action_base, mode, lambda_te, s, t )

if (t > 1)
    % we only select actions, after the first iteration, so that we
    % first define the baseline action and utility
    if ( mode(s,1) == 1 )
        % if this player is 'content'
        P_aggre = zeros(No_actions(s,1),1);
        sum_p = 0;  % sum of probability entries
        for j = 1 : No_actions(s,1)
            if ( j ~= action_base(s,1) )
                % with probability epsilon, it plays any other 
                % action with equal probability each 
                P_aggre(j,1) = lambda_te/(No_actions(s,1)-1) + sum_p;
                sum_p = P_aggre(j,1);		% update of the sum
            else
                % with probability 1 - epsilon, it plays the same
                % action
                P_aggre(j,1) = 1 - lambda_te + sum_p;
                sum_p = P_aggre(j,1);
            end
        end
        action_out = mixed_strategy(P_aggre',No_actions(s,1));
    elseif ( mode(s,1) == 0 )
        % if the player is discontent, then it plays any action
        % uniformly at random.
        P_aggre = zeros(No_actions(s,1),1);
        sum_p = 0;  % sum of probability entries
        for j = 1 : No_actions(s,1)
            P_aggre(j,1) = 1/No_actions(s,1) + sum_p;
            sum_p = P_aggre(j,1);		% update of the sum
        end
        action_out = mixed_strategy(P_aggre',No_actions(s,1));
    elseif ( mode(s,1) == 2 )
        % if the player is watchful, then it plays the benchmark
        % action
        action_out = action_base(s,1);
    else
        % if player is hopeful
        action_out = action_base(s,1);
    end
else
    action_out = action_base(s,1);
end


end

