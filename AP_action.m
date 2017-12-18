function action_out = AP_action ( action_prev, No_actions, No_neig, Neighbors, con_symaloc, DIAM, ns, s, ...
    rho, decay, kappa_0, kappa_1, distances_global, No_down_nodes, down_nodes, deficiency, const_def, P, R, BR, Neigh_payoff, p_ada)

% na            = No_actions(s,1)
% 

R_trial = zeros(No_actions(s,1),1); % vector of rewards over actions
for j=1:No_actions(s,1) % for each available action
    % setting up trial action
    action_trial = action_prev;
    action_trial(s,1)=j;
    %--------------------
    % Computing Rewards
    % reward preliminaries
    [No_links,L,distances,distances_global,down_nodes,No_down_nodes,LtN] = ...
        reward_preliminaries(ns,DIAM,No_neig,Neighbors,No_actions,action_trial);
    % reward
    R = reward(ns,DIAM,rho,decay,L,kappa_0,kappa_1,No_neig,Neighbors,con_symaloc,...,
        distances_global,No_down_nodes,down_nodes,deficiency,const_def,No_links,...
        No_actions,action_trial,P,Neigh_payoff);
    R_trial(j,1) = R{s,1};
end
% computing better reply
BR{s,1}=0;  % better reply initialization
[R_trial_max,BR{s,1}(1,1)] = max(R_trial); % sorted in ascending order
ccc=1; % counter
for j = 1:No_actions(s,1)
    if R_trial(j,1) == R_trial_max & j ~= BR{s,1}(1,1)
        % if the value equal the maximum value
        ccc=ccc+1;
        BR{s,1}(ccc,1) = j;
    end
end
% checking if better-reply is empty
BR_empty = 0;
for j = 1:ccc
    if BR{s,1}(j,1) == action_prev(s,1)
        BR_empty = 1;   % better reply is empty
    end
end

% selecting action
if BR_empty==1 
    action(s,1) = action_prev(s,1);
elseif BR_empty==0
    temp=rand(1);
    if temp < p_ada
        action(s,1) = action_prev(s,1);
    else
        action(s,1) = random_selection(BR{s,1}); % randomly select one of the better replies
    end
end

action_out = action(s,1);

end

