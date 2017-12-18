function [mode_new,action_base_new,R_base_new] = TE_update(R,R_base,action,action_base,mode)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    if ( mode == 1 )
        % if currently the player is content
        if ( (action == action_base) && (R == R_base) )
            % if its action and payoff remained the same, it remains
            % content, the base action/payoff do not change
            mode_new = 1;
            action_base_new = action_base;
            R_base_new = R_base;
        elseif ( action ~= action_base && R <= R_base )
            % if the action changed and the payoff is less or equal to
            % the base payoff, then nothing changes
            mode_new = 1;
            action_base_new = action_base;
            R_base_new = R_base;
        elseif (action ~= action_base && R > R_base )
            mode_new = 1;
            action_base_new = action;
            R_base_new = R;
        elseif (action == action_base && R < R_base )
            mode_new = 2;  % watchful
            action_base_new = action_base;
            R_base_new = R_base;
            % benchmark action/payoff do not change
        elseif (action == action_base && R > R_base )
            mode_new = 3;  % hopeful
            action_base_new = action_base;
            R_base_new = R_base;
        end
    elseif (mode == 2)
        % watchful
        if (action == action_base && R < R_base)
            mode_new = 0;
        elseif (action == action_base && R == R_base)
            mode_new = 1;
        elseif (action == action_base && R > R_base) 
            mode_new = 3;
        end
        action_base_new = action_base;
        R_base_new = R_base;
    elseif (mode == 3)
        % hopeful
        if (action == action_base && R < R_base)
            mode_new = 2;
            R_base_new = R_base;
        elseif (action == action_base && R == R_base)
            mode_new = 1;
            R_base_new = R_base;
        elseif (action == action_base && R > R_base)
            mode_new = 1;
            R_base_new = R;
        end
        action_base_new = action_base;
    elseif (mode == 0)
        % discontent
        P_aggre = zeros(2,1);
        % probability of becoming discontent
        P_aggre(1,1) = 1 - sigmoid_Young(R - R_base,0.01);
        P_aggre(2,1) = 1;
        mode_new = mixed_strategy(P_aggre',2) - 1;
        % randomly choose between content/discontent
        if (mode_new == 1)
            % if content, then we update the benchmarks based on the
            % current action/payoff
            action_base_new = action;
            R_base_new = R;
        elseif (mode_new == 0)
            action_base_new = action_base;
            R_base_new = R_base;
        end
    end


end

