function [mode_new,action_base_new,R_base_new] = MB_update(R,R_base,action,action_base,mode,lambda_mb)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

    if ( mode == 1 )
        if ( (action == action_base) && (R == R_base) )
            mode_new = 1;
            R_base_new = R_base;
            action_base_new = action_base;
        else
            action_base_new = action;
            R_base_new = R;
            % Randomly selecting mode
            P_aggre = zeros(2,1);
            P_aggre(1,1) = 1 - lambda_mb^(1-R);
            P_aggre(2,1) = 1;
            mode_new = mixed_strategy(P_aggre',2) - 1;
        end


    else % if it is discontent

        action_base_new = action;
        R_base_new = R;
        % Randomly selecting mode
        P_aggre = zeros(2,1);
        P_aggre(1,1) = 1 - lambda_mb^(1-R);
        P_aggre(2,1) = 1;
        mode_new = mixed_strategy(P_aggre',2) - 1;
    end
        
end

