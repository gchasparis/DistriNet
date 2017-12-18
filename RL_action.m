function [ P, action ] = RL_action ( epsilon, P, P_nor, P_sel, Q, g, na, PROJ )

    % na        = No_actions(s,1)
    % P         = P{s,1}
    % Q         = Q{s,1}
    % epsilon   = epsilon
    % P_nor     = P_nor{s,1}
    % g         = g(gain,1)
    
    %-------------------
    % DERIVATIVE ACTION
    P = P + g * (P - Q);
    
    %-------------------
    % PROJECTION
    if PROJ==1
        
        %     sum_p = sum(P{s,1});
        %     c = 0; %counter of negative elements
        %     for i = 1 : No_actions(s,1)
        %         if P{s,1}(i,1) < 0;
        %             c = c+1;
        %             break
        %         end
        %     end
        % 
        %     if  abs(sum_p-1) > 10e-6 | c >= 1   
        %         P{s,1} = projection_to_simplex_lsq(P{s,1},No_actions(s,1));
        %     end

        % Alternative (suboptimal projection)
        for i = 1 : na
            if P(i,1) < 0;
                P(i,1) = 0;
            elseif P(i,1) > 1;
                P(i,1) = 1;
            end    
        end

    end % of PROJ
    
    %-----------------
    % PICK AN ACTION
    P_sel = (1 - epsilon)* P + epsilon * P_nor;

    % Create aggregate probability vector
    P_aggre = zeros(na,1);
    sum_p = 0;  % sum of probability entries
    for j = 1 : na
        P_aggre(j,1) = P_sel(j,1) + sum_p;
        sum_p = P_aggre(j,1);		% update of the sum
    end
    
    % Players first choose an action
    action = mixed_strategy(P_aggre',na);
    
end

