function [No_links,L,distances,distances_global,down_nodes,No_down_nodes,LtN] = ...
    NF_reward_preliminaries(ns,DIAM,No_neig,Neighbors,No_actions,action)


LtN=cell(ns,1);
distances=cell(ns,1);
distances_global=cell(ns,1);
down_nodes=cell(ns,1);
No_down_nodes=cell(ns,1);

%----------------------
% COMPUTATION OF LINKS
%----------------------
for s=1:ns
    AtL = actions_to_links(No_neig(s,1),No_actions(s,1),1);
    LtN{s,1} = links_to_neighbors(AtL,Neighbors{s,1},action(s,1),s);
    if LtN{s,1}(1,1) ~= s    % if you don't connect with noone
        [No_links(s,1),temp] = size(LtN{s,1});
    else
        No_links(s,1) = 0;
    end
        
end % of s

%-------------------------
% COMPUTATION OF LAPLACIAN
%-------------------------
% Laplacian Matrix (source of information - from row to column)
L = zeros(ns,ns);
for s = 1 : ns
    for n = 1 : No_links(s,1)
        L(LtN{s,1}(n,1),s) = 1;
    end
end


%--------------------------
% COMPUTATION OF DISTANCES 
%--------------------------
% computation of distances from neighbors
for s = 1:ns
    distances{s,1}=zeros(No_neig(s,1)-1,1);%inf*ones(No_neig(s,1)-1,1);
%     connectivity{s,1}=zeros(No_neig(s,1)-1,1);
    temp=zeros(No_neig(s,1)-1,1);%inf*ones(No_neig(s,1)-1,1);
    for j = 1 : No_neig(s,1)-1
        % for each neighbor
        counter=0;
        for m = 1 : ns-1
            XX=(L')^m;
            if XX(s,Neighbors{s,1}(j,1))>0 & counter==0
                counter=counter+1;
                temp(j,1)=m;
            end
        end
        distances{s,1}(j,1) = temp(j,1);
    end
end

% computation of distances from everybody
for s = 1:ns
    distances_global{s,1}=zeros(ns,1);
    temp=zeros(ns,1);
    for j=1:ns
        if j ~= s   % for every other node
            counter=0;
            for m = 1 : ns-1
            XX=(L')^m;
            if XX(s,j)>0 & counter==0
                counter=counter+1;
                temp(j,1)=m;
            end
            end
        end
        distances_global{s,1}(j,1) = temp(j,1);
    end
end
            

%---------------------------------
% COMPUTATION OF DOWNSTREAM NODES
%---------------------------------
if DIAM>0   % if there is a constraint on the diameter
    for s=1:ns  % for each node
        down_nodes{s,1}=zeros(1,1);
        No_down_nodes{s,1}=0;
        
        for j=1:ns
            if j ~= s & distances_global{j,1}(s,1) <= DIAM-1 & ...
                    distances_global{j,1}(s,1) > 0 % for each other node within distance DIAM
                No_down_nodes{s,1} = No_down_nodes{s,1}+1;
                down_nodes{s,1}(No_down_nodes{s,1},1) = j;
            end
        end        
    end
end
