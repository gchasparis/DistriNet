function [ ASP_out ] = AL_update ( asp_pert, ns, ASP, R, s, zeta, asp_min, asp_max, com_pert, mu_asp )

ASP_out = ASP;

if com_pert==0

    for s=1:ns
        xx=rand(1);
        if xx> asp_pert
            ASP_out{s,1} = ASP{s,1} + mu_asp*(R{s,1}-ASP{s,1});
            ASP_out{s,1} = sat(ASP{s,1},asp_min,asp_max);
        else
            % perturbation = random('Normal',0,zeta,1,1);
            % ASP{s,1} = ASP{s,1} + mu_asp*(R{s,1}-ASP{s,1}) + perturbation; 2*zeta*rand(1,1) - zeta
            ASP_out{s,1} = ASP{s,1} + mu_asp * (R{s,1}-ASP{s,1}) + 2 * zeta * rand(1,1) - zeta;
            ASP_out{s,1} = sat(ASP{s,1},asp_min,asp_max);
        end
    end

else % players use combined perturbation

    xx=rand(ns,1);
    theta=zeros(ns,1);
    for s = 1:ns
        if xx(s) > asp_pert
            theta(s)=0;
        else
            theta(s)=1;
        end
    end
    % compute maximum perturbation
    theta_max = max(theta(:,1));

    % computation of perturbation bounds
    xxx_up=zeros(ns,1); % upper bound on perturbation
    xxx_down=zeros(ns,1); % lower bound on perturbation
    for s=1:ns
        % lower perturbation bound
        if ASP{s,1}-asp_min <= zeta
            xxx_down(s,1) = ASP{s,1}-asp_min;
        else
            xxx_down(s,1) = zeta;
        end
        % upper perturbation bound
        if asp_max-ASP{s,1} <= zeta
            xxx_up(s,1) = asp_max-ASP{s,1};
        else
            xxx_up(s,1) = zeta;
        end
    end

    % actual perturbation
    xxx=zeros(ns,1); % actual perturbation
    for s=1:ns
        if theta_max==0
            xxx(s,1)=0;
        elseif theta_max==1 & theta(s)==0
            xxx(s,1)=unifrnd(-xxx_down(s,1),0);
        elseif theta(s)==1
            xxx(s,1)=unifrnd(-xxx_down(s,1),xxx_up(s,1));
        end
        % aspiration update
        ASP_out{s,1} = ASP{s,1} + mu * (R{s,1}-ASP{s,1}) + xxx(s,1);
    end

end % of combined pertubation check

end

