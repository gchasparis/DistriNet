function LtN = links_to_neighbors(AtL,Neighbors,action,s);

No_links = AtL(action,1);

LtN = zeros(No_links,1);  % number of links

for i = 1 : No_links
    temp = AtL(action,i+1);
    if temp ~= 1        % you have a link with yourself
        LtN(i,1) = Neighbors(AtL(action,i+1)-1,1);
    else
        LtN(i,1) = s;   % yourself
    end        
end
    