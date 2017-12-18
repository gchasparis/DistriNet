function R = NF_reward(ns,DIAM,rho,decay,L,kappa_0,kappa_1,No_neig,Neighbors,con_symaloc,...
    distances_global,No_down_nodes,down_nodes,deficiency,const_def,No_links,No_actions,action,P,Neigh_payoff);

R = cell(ns,1);

temp1 = cell(ns-1,ns);
temp1_nodecay = cell(ns-1,ns);
mai_cost = zeros(ns,1);
est_cost = zeros(ns,1);

%---------------    
if DIAM==0 & Neigh_payoff==0 % if there is no constraints on the DIAMETER and we only consider Directional Links
    
    for s = 1 : ns  % for all nodes

    
%     % it can be computed based only on distances
%     
%     R_test{s,1}=1;   % initialization (benefit from itself)
%     
%     for j = 1 : No_neig(s,1)-1  % from all s's neighbors
%         if distances{s,1}(j,1) < Inf
%             R_test{s,1} = R_test{s,1} + 1;
%         end
%     
    
    %--------
    R{s,1} = rho * decay * ones(1,ns) * L * vertex(ns,s);
    
    for h = 1 : ns-1; % maximum distance between any two nodes
        temp1{h,s} = decay^h * L^h * vertex(ns,s);
        temp1_nodecay{h,s} = L^h * vertex(ns,s);
        %-------------
        % make sure you don't count a node twice
        for i = 1 : ns
            if temp1_nodecay{h,s}(i) ~= 0 
                temp1_nodecay{h,s}(i) = 1;
                temp1{h,s}(i) = decay^h;
            end
        end   
    end
    %-------------
    
    for j = 1 : (ns-1)-1
        temp2 = vertex(ns,s);   % by defining it this way, we exclude the 
        % information stemming from the current node 
        for k = 1 : j
            temp2 = temp2 + temp1_nodecay{k,s};
        end
        for n = 1 : ns
            if temp2(n) ~= 0
                temp2(n) = 1;
            end
        end

        R{s,1} = R{s,1} + rho * (ones(ns,1) - temp2)' * temp1{j+1,s};
    end

    %------------
    % cost of maintaining a link
    mai_cost(s,1) = kappa_0 * No_links(s,1);

    % cost of establishing new link
    a{s,1} = vertex(No_actions(s,1),action(s,1));
    est_cost(s,1) = kappa_1 * No_links(s,1) * (ones(1,No_actions(s,1))-a{s,1}')...
        * P{s,1};
    %------------
    
    % total benefits
    R{s,1} = R{s,1} - est_cost(s,1) - mai_cost(s,1);

    end % of s
    
   

%--------------
elseif DIAM==0 & Neigh_payoff==1
 
    for s = 1 : ns  % for all nodes

    R{s,1}=0; %initialization
    for j=1:No_neig(s,1)-1
        % for all neighbors of agent s
        if distances_global{s,1}(Neighbors{s,1}(j,1),1) > 0
            R{s,1} = R{s,1} + 1;
        end
    end

    %------------
    % cost of maintaining a link
    mai_cost(s,1) = kappa_0 * No_links(s,1);

    % cost of establishing new link
    a{s,1} = vertex(No_actions(s,1),action(s,1));
    est_cost(s,1) = kappa_1 * No_links(s,1) * (ones(1,No_actions(s,1))-a{s,1}')* P{s,1};
    %------------
    
    % total benefits
    R{s,1} = R{s,1} - est_cost(s,1) - mai_cost(s,1);

    end % of s
    
    %disp '-------------------'
    %R
    
    total_R=0;
    for s=1:ns
        total_R = total_R + R{s,1};
    end
        
    if con_symaloc == 1
        for s=1:ns
            R{s,1} = 1/ns * total_R;
        end
        %R
    end

    
%--------------    
elseif DIAM~=0 % if there is a constraint on the DIAMETER
    
    % first we compute production/deficiency for each node
    for s = 1 : ns
        
    %-------
    % first we computer the production for each node
    % (we define the production of agent s as the benefits it 
    % receives from its neighbors within distance DIAM)
    production{s,1}=0;
    for j = 1 : No_neig(s,1)-1 % for each neighbor
        if distances_global{s,1}(Neighbors{s,1}(j,1),1) <= DIAM & ...
                distances_global{s,1}(Neighbors{s,1}(j,1),1) > 0
            production{s,1} = production{s,1} + 1;
        else
            production{s,1} = production{s,1};
        end
    end
    % Deficiency of downstream (within DIAM) nodes
    % (we define the deficiency of agent s as the number of neighbors
    % that are not being accessed)
    deficiency{s,1}=0; %if there is no deficiency (all neighbors are accessed)
    for j = 1 : No_neig(s,1)-1 % for each neighbor of s
        if distances_global{s,1}(Neighbors{s,1}(j,1),1) > DIAM | ...
                distances_global{s,1}(Neighbors{s,1}(j,1),1) == 0
            % if the distance from j to s is greater than DIAM (very far) or 0 
            % (there is no link)
            % (i.e., if s is not accessing j)
            deficiency{s,1} = deficiency{s,1} + 1;
        end
    end
    
    % Number of downstream deficient nodes
    No_down_deficients{s,1} = 0; % initialization
    if No_down_nodes{s,1}>0 % if there are downstream nodes
        for j = 1 : No_down_nodes{s,1};
            % for each of the nodes downstream
            if deficiency{down_nodes{s,1}(j,1),1}>0
                No_down_deficients{s,1} = No_down_deficients{s,1} + 1;
            end
        end
    end
            
    
    %-------
    % reward based on downstream deficients
    temp = production{s,1} * (1 - const_def * No_down_deficients{s,1}/(1+No_down_deficients{s,1}));
    
    %------------
    % cost of maintaining a link
    mai_cost(s,1) = kappa_0 * No_links(s,1);

    % cost of establishing new link
    a{s,1} = vertex(No_actions(s,1),action(s,1));
    est_cost(s,1) = kappa_1 * No_links(s,1) * (ones(1,No_actions(s,1))-a{s,1}') * P{s,1};
    %------------

    R{s,1} = temp - est_cost(s,1) - mai_cost(s,1);
    
    
%     %-------
%     % reward of agent s based on downstream productions
%     % (it is defined as the avergage production of all downstream nodes)
%     temp = production{s,1};
%     % reward based on downstream productions
%     if No_down_nodes{s,1}>0 % if there are downstream nodes
%     for j = 1 : No_down_nodes{s,1} % for all nodes downstream
%         temp = temp + production{down_nodes{s,1}(j,1),1};
%     end
%     end
%     
%     %------------
%     % cost of maintaining a link
%     mai_cost(s,1) = kappa_0 * No_links(s,1);
% 
%     % cost of establishing new link
%     a{s,1} = vertex(No_actions(s,1),action(s,1));
%     est_cost(s,1) = kappa_1 * No_links(s,1) * (ones(1,No_actions(s,1))-a{s,1}')...
%         * P{s,1};
%     %------------
%     
%     R{s,1} = 1/(No_down_nodes{s,1}+1) * temp - est_cost(s,1) - mai_cost(s,1);
    
    
% elseif DIAM~=0 & Neigh_payoff==1
%     % if there is a constraint on the diameter and the payoff is 
%     % based only on the neighborhood
% 
%     % In this case, the payoff will coincide with the production
%     
%     %-------
%     % first we computer the production for each node
%     % (we define the production of agent s as the benefits it 
%     % receives from its neighbors within distance DIAM
%     production{s,1}=0;
%     for j = 1 : No_neig(s,1)-1 % for each neighbor
%         if distances{s,1}(j,1) <= DIAM
%             production{s,1} = production{s,1} + 1;
%         else
%             production{s,1} = production{s,1};
%         end
%     end
%     
%     %------------
%     % cost of maintaining a link
%     mai_cost(s,1) = kappa_0 * No_links(s,1);
% 
%     % cost of establishing new link
%     a{s,1} = vertex(No_actions(s,1),action(s,1));
%     est_cost(s,1) = kappa_1 * No_links(s,1) * (ones(1,No_actions(s,1))-a{s,1}')...
%         * P{s,1};
%     %------------
% 
%     R{s,1} = production{s,1} - est_cost(s,1) - mai_cost(s,1);


    end % of s


end % of if DIAM

% normalizing reward [0,1]
for (s = 1 : ns)
    R{s,1} = R{s,1} / ( ns - 1 ) / rho; 
end
