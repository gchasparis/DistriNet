function AtL = actions_to_links(np,na,row)

%np=4;na=8;row=1;

% we introduce a vector of players
vector_players = zeros(np,1);
for i = 1 : np
    vector_players(i,1) = i;
end

% we introduce a matrix that assigns actions to links
AtL = zeros(na,np-1+1);

AtL(1:np,1:2) = [ones(np,1) vector_players];
n_comb_tot = np;     % counter
for k = 2 : np-1
    % here we remove the player row from the set
    vector_players_reduced = zeros(np-1,1);
    c = 0;  % counter
    for h = 1 : np
        if h ~= row
            c = c + 1;
            vector_players_reduced(c,1) = vector_players(h,1);
        end
    end
    comb = nchoosek(vector_players_reduced,k);
    [n_comb,x] = size(comb);
    AtL( (n_comb_tot+1) : (n_comb_tot+n_comb),1:k+1) = [k*ones(n_comb,1) comb];
    n_comb_tot = n_comb_tot + n_comb;
end


