% /* -----------------------------------------------------------------------------
%  * Copyright (C) 2012-2014 Software Competence Center Hagenberg GmbH (SCCH)
%  * <georgios.chasparis@scch.at>, <office@scch.at>
%  * -----------------------------------------------------------------------------
%  * This program is free software: you can redistribute it and/or modify
%  * it under the terms of the GNU General Public License as published by
%  * the Free Software Foundation, either version 3 of the License, or
%  * (at your option) any later version. 
%  *
%  * This program is distributed in the hope that it will be useful,
%  * but WITHOUT ANY WARRANTY; without even the implied warranty of
%  * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  * GNU General Public License for more details.
%  *
%  * You should have received a copy of the GNU General Public License
%  * along with this program.  If not, see <http://www.gnu.org/licenses/>.
%  * -----------------------------------------------------------------------------
%  * Upon use of the provided software, it is kindly requested that the 
%  * following scientific publication is cited (based on which this 
%  * software has been developed):
%  * 
%  * Georgios Chasparis, Jeff Shamma, "Network Formation: Neighbohood
%  * Structures, Establishment Costs, and Distributed Learning", IEEE
%  * Transactions on Cybernetics, vol 43, no 6, 2013.
%  * (http://ieeexplore.ieee.org/document/6425448/)
%  * 
% */

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description

% This software simulation addresses distributed learning in a network
% formation game. This simulator is based upon the network formation model
% presented in reference [1] (see details below).
% 
% In this software, we allow the user to select alternative payoff-based
% learning dynamics. In particular, we have already incorporated the
% following classes of payoff-based learning:
% 
% [RL]: Reinforcement Learning (see reference [1])
% [AR]: Apiration-based Reinforcement Learning (see references [2], [3])
% [AL]: Aspiration Learning (see references [4], [5])
% [AP]: Adaptive Play (see reference [6])
% [BB]: Baseline-based Learning (see reference [7])
% [MB]: Mode-based Learning (see reference [8])
% [TE]: Trial-and-Error Learning (see reference [9])
% 
% References:
%
% [1] Georgios Chasparis, Jeff Shamma, "Network Formation: Neighbohood
% Structures, Establishment Costs, and Distributed Learning", IEEE
% Transactions on Cybernetics, vol 43, no 6, 2013.
% [2] Georgios Chasparis, Michael Rossbory, Vladimir Janjic, "Efficient
% Dynamic Pinning of Parallelized Applications with Reinforcement Learning
% with Applications," Lecture Notes in Computer Science, Euro-Par 2017:
% Parallel Processing, vol 10417, 2017.
% [3] Georgios Chasparis, "Aspiration-based Learning Automata", (submitted
% to) European Control Conference 2018.
% [4] Georgios Chasparis, Ari Arapostathis and Jeff Shamma, "Aspiration
% Learning in Coordination Games", SIAM Journal on Control and
% Optimization, vol 51, no 1, 2013.
% [5] R. Karandikar, D. Mookherjee, and D. Ray, "Evolving Aspirations and
% Cooperation," Journal on Economic Theory, vol 80, 1998.
% [6] P. Young, "The evolution of conventions," Econometrica, vol 61, pp
% 57-84 1993.
% [7] J. Marden, H. P. Young, G. Arslan and J. Shamma, "Payoff based 
% dynamics for multi-player weakly acyclic games," SIAM Journal on Control 
% and Optimization, vol 48, no 1, pp. 373--396.
% [8] J. Marden, H. P. Young and L. Pao, "Achieving Pareto Optimality
% through Distributed Learning," vol 52, no 5, pp. 2753-2770, 2014.
% [9] P. Young, "Learning by trial and error", Games and Economic Behavior,
% vol 65, no 2, pp 626-643, 2009.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;

%% Setup
% Payoff-based Learning Dynamics 
% Please set 1 for the preferred type of learning dynamics
% (AR = 1 should be combined with RL = 1)
RL = 1;                                 % Reinforcement Learning 
AR = 1;                                 % Aspiration-based Reinforcement Learning (applies if and only if RL = 1 as well)
AL = 0;                                 % Aspiration Learning (applies by setting AL = 1)
AP = 0;                                 % Adaptive Play (applies by setting AP = 1)
BB = 0;                                 % Baseline-based Learning (applies by setting BB = 1)
MB = 0;                                 % Mood-based Learning (applies by setting MB = 1)
TE = 0;                                 % Trial-and-Error Dynamics (applies by setting TE = 1)

% Visualization Tools 
graph_plot      = 0;                    % plotting the graph 
time_plot       = 0;                    % plotting responses with time
aver_plot       = 0;                    % plotting average performances with time
dist_plot       = 0;                    % if you want to plot distances set "1" (otherwise "0")
camera          = 0;                    % recording network evolution
strategies_plot = 0;                    % plotting strategy vectors with time (applies to RL and AR)
mu_plot         = 0.1;                  % step-size for running average behavior
size_noise      = 0.02;                 % size of noise term in observations

% Network Formation Reward Parameters (check also reference [1] above)
loc = 1;                                % if you only count distances from neighbors (otherwise set "0")
DIAM = 0;                               % diameter constraint. "0" if none
Neigh_payoff = 0;                       % payoff dependent only on neighbors
con_symaloc = 0;                        % if the payoff is based on a symmetric allocation of the connections model
                                        % need to also set DIAM=0.
rho             = 4;                    % reward per link
kappa_0         = 2;                    % maintenance cost
kappa_1         = 0;                    % establishment cost
decay           = 1;                    % decay factor
const_def       = 1/2;                  % constant for reward based on deficients
                        
% Generic Parameters (e.g., storage)
No_runs = 10;                           % (repeat execution for multiple times)
flash_interval = 20;                    % frequency of capturing frame
rect = [72 44 440 346];                 % for camera frame
T = 5000;                               % number of loops per simulation run
t_stor_per = 50;                        % frequency for storing variables
T_stor = (T/t_stor_per);
T_stor_vec(1,1) = 0;


%% Payoff-based Learning Parameters

% Reinforcement Learning (RL) Parameters (check reference [1])
aver            = 0;                        % for decreasing step-size select "aver=1"
q               = 7/9;                      % decreasing step-size exponent (1/(t^q+1))
mu              = 0.2;                      % step size for reinforcement learning when "aver=0" (constant step-size)
epsilon         = 0.01;                     % mutation rate (probability with which a uniform distribution is applied in the action selection)
T_aver          = 100;                      % iteration after which a decreasing step-size is implemented (if such selected)
p_rl            = 0.01;                     % sigmoid parameter under AR=1
c_rl            = 1/p_rl;                   % sigmoid parameter under AR=1

% Dynamic Reinforcement Update (check reference [1])
No_gains = 1;
g = [0.01];                                 % example [0.0;1/200;1/150;1/100;1/50];
PROJ=0;                                     % Projection on the simplex (usually required under Dynamic Reinforcement)

% Aspiration Learning (AL) Parameters (check reference [4])
p_asp           = 0.01;                     % aspiration minimum probability
c_asp           = 5;                        % this is the slope of the sigmoid function
asp_pert        = 0.01;                     % probability of perturbation of the aspiration level update (lambda)
mu_asp          = 0.02;                     % step-size of aspiration learning
zeta            = 0.01;                     % this is the size of noise in the aspiration-level update
com_pert        = 0;                        % if 1 combined perturbation is applied
asp_max         = 1.01;                     % upper bound of aspiration levels
asp_min         = -0.01;                    % lower bound of aspiration levels
asp_sta         = 0.0;                      % starting aspiration level

% Adaptive Play (AP) Parameters 
M = 2;                  % memory size
p_ada = 0.1;

% Baseline-Based Learning (BB) Parameters 
epsilon_bb = 0.01;

% Mood-based Learning (MB) Parameters 
lambda_mb = 0.1;
c_mb = 17;

% Trial-and-Error Learning (TE) Parameters 
lambda_te = 0.5;


%% Grid Parameters & Setup
% Grid Details 
ns = 16;                                    % number of sensors 16 for regular grid
r_in = .01;                                 % radius of sensor
r_out = .31;                                % radius of sensing area 0.7 (.21 when you use 17 nodes test_grid)
                                            % 0.38 for test_grid_3
                                            % 0.6 for the simple test of 6 nodes
arrow_tipangle  = 10;
arrow_length    = 14;

load grid/test_grid x y                     % loading nodes coordinates

% Creating Graph
x_c_in = cell(ns,1);
y_c_in = cell(ns,1);
x_c_out = cell(ns,1);
y_c_out = cell(ns,1);
for s = 1 : ns
k = 0;      % counter
for theta = 0 : 0.01 : 2*pi
    k = k + 1;
    % object limits
    x_c_in{s,1}(k,1) = x(s,1) + r_in * cos(theta);
    y_c_in{s,1}(k,1) = y(s,1) + r_in * sin(theta);
    
    % vision limits
    x_c_out{s,1}(k,1) = x(s,1) + r_out*cos(theta);
    y_c_out{s,1}(k,1) = y(s,1) + r_out*sin(theta);
end
end

%Setting up neighbors
Neighbors = cell(ns,1);
for s = 1 : ns  % for any sensor
    Neighbors{s,1} = 0;
    k = 0;
    for b = 1 : ns % for any other sensor
        if s ~= b & sqrt((x(s,1)-x(b,1))^2 + (y(s,1)-y(b,1))^2)<r_out
            k = k + 1;
            Neighbors{s,1}(k,1) = b;
        end
    end
end
            
No_neig = zeros(ns,1);
for s = 1 : ns
    if Neighbors{s,1}(1,1) ~= 0
        [No_neig(s,1),temp] = size(Neighbors{s,1});   
        No_neig(s,1) = No_neig(s,1)+1;   % including the sensor itself
    elseif Neighbors{s,1}(1,1) == 0
        No_neig(s,1) = 1;
    end
end

No_actions = zeros(ns,1);
for s = 1 : ns
    No_actions(s,1) = No_neig(s,1);
    for k = 2 : No_neig(s,1)-1
        No_actions(s,1) = No_actions(s,1) + factorial(No_neig(s,1)-1)/...
            (factorial(k)*factorial(No_neig(s,1)-1-k));
    end
end


%% Line Formats (in the output figures)
S = ['k- ';'k--';'k-.';'kv ';'kd ';'k* ';'ko ';'k+ ';'kx ';'k^ ';'kp ';'ks ';...
        'k. ';'k< ';'k> ';'kh ';'k* '];
S_col = ['b  ';'r  ';'g  ';'k  ';'b--';'r--';'g--';'k--';'b: ';...
   'r: ';'g: ';'k: ';'b-.';'r-.';'g-.';'k-.';'b. ';'r. ';'g. ';'k. '];

% Legends
A = ['Action A';'Action B';'Action C';'Action D';'Action E';...
    'Action F';'Action G';'Action H';'Action I';'Action J';...
    'Action K';'Action L';'Action M';'Action N';'Action O';...
    'Action P'];

Agent = ['Agent  1';'Agent  2';'Agent  3';'Agent  4';'Agent  5';...
    'Agent  6';'Agent  7';'Agent  8';'Agent  9';'Agent 10';...
    'Agent 11';'Agent 12';'Agent 13';'Agent 14';'Agent 15';...
    'Agent 16';'Agent 17';'Agent 18';'Agent 19';'Agent 20'];

G = ['\gamma_{1}';'\gamma_{2}';'\gamma_{3}';'\gamma_{4}';'\gamma_{5}'];

Trial = ['./Movies/Trial_0';
    './Movies/Trial_1';
    './Movies/Trial_2';
    './Movies/Trial_3';
    './Movies/Trial_4';
    './Movies/Trial_5';
    './Movies/Trial_6';
    './Movies/Trial_7';
    './Movies/Trial_8';
    './Movies/Trial_9'];

gain_legend_array = cellstr(G);
agent_legend_array = cellstr(Agent);
legend_array = cellstr(A);
line_format = cellstr(S);
line_format_col = cellstr(S_col);

Trials_array = cellstr(Trial);


%% Running-Average Storage Parameters 
% (Initialization of running average performance counters)
R_ave = cell(No_gains,1);
D_ave = cell(No_gains,1);
R_run_ave = cell(No_gains,No_runs);
R_run_ave_sto = cell(No_gains,No_runs);
R_run_ave_per_run = cell(No_gains,1);


%% Loop (for Dynamic Reinforcement Gain if used)
for gain = 1 : No_gains

R_ave{gain,1} = zeros(No_runs,1);
D_ave{gain,1} = zeros(No_runs,1);
R_run_ave_per_run{gain,1} = zeros(No_runs,1);

%% Loop (for each Run)
for run = 1 : No_runs
    
disp '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~';
[gain,run]

if camera == 1;
    avio = avifile(char(Trials_array(run)),'compression','Cinepak','fps',1);
end

% Definition of storage variables of Payoff-based Learning Dynamics
if RL==1
    % Definitions on Reinforcement Learning
    P = cell(ns,1);     % Probability distributions of actions for all players
    Q = cell(ns,1);     % Previous probability distribution
    Y = cell(ns,1);     % A new state that tracks rewards
    r = cell(ns,1);     % Vectors of the rewards
    R = cell(ns,1);     % Aggregate rewards
    R_sto = cell(ns,1); % Storing rewards
    action_vec = cell(ns,1);

    P_nor = cell(ns,1); % Normal probability distribution
    P_sel = cell(ns,1); % Probability distribution for selection
    P_sto = cell(ns,1); % stored probability distribution
    Y_sto = cell(ns,1); % stored value of Y state
elseif AL==1
    % Definitions on Aspiration Learning
    ASP = cell(ns,1);   % Aspiration level of all agents
    ASP_sto = cell(ns,1);% Stored aspiration level of all agents
    P = cell(ns,1);     % probability of selecting actions
elseif AP==1
    %Definitions on Adaptive Play
    P = cell(ns,1); % ACTION = Empirical Frequency of Actions over Finite Interval M
    action_vec = cell(ns,1); % ACTION = The current action 
    action_vec_hist = cell(ns,1); % History of Actions over finite interval M
    % Storing variables
    P_sto = cell(ns,1); % Storing the state vector
elseif BB==1
    % Definitions on Benchmark-Based Dynamics
    R = cell(ns,1);             % Aggregate rewards
    R_base = cell(ns,1);
    R_sto = cell(ns,1);         % Storing rewards
    action_vec = cell(ns,1);
    action_vec_hist = cell(ns,1);
    action_base = zeros(ns,1);
    P = cell(ns,1);     % probability of selecting actions (not needed)
elseif MB==1
    R = cell(ns,1);
    R_base = cell(ns,1);
    R_sto = cell(ns,1);
    action_vec = cell(ns,1);
    action_vec_hist = cell(ns,1);
    action_base = zeros(ns,1);
    mode = zeros(ns,1);
    mode_base = zeros(ns,1);    % 0 for discontent and 1 for content
    P = cell(ns,1);     % probability of selecting actions (not needed)
elseif TE==1
    R = cell(ns,1);
    R_base = cell(ns,1);
    R_sto = cell(ns,1);
    action_vec = cell(ns,1);
    action_vec_hist = cell(ns,1);
    action_base = zeros(ns,1);
    mode = zeros(ns,1);         % 0: discontent, 1: content, 2: watchful, 3: hopeful
    P = cell(ns,1);
end

if dist_plot==1
    % Distances from neighbors (definition)
    distances = cell(ns,1);
    distances_sto = cell(ns,1);
    distances_ave = cell(ns,1); % average distances over time
    distances_ave_sto = cell(ns,1);
    distances_inf = cell(ns,1);
    inv_distances_run_ave = zeros(ns,1);
    inv_distances_run_ave_sto = cell(ns,1);
    total_inv_distances_run_ave = 0;
    total_inv_distances_run_ave_sto = zeros(T_stor,1);
end
% Distances from everybody (definition)
distances_global = cell(ns,1);

% Production
if DIAM>0;
    production = cell(ns,1);
    deficiency = cell(ns,1);
    No_down_deficients = cell(ns,1);
    down_nodes = cell(ns,1);
    No_down_nodes = cell(ns,1);
end


% General Definitions
LtN = cell(ns,1);   % Neighbors that correspond to the current links

%----------------------------------------------
% Initializing parameters
action = zeros(ns,1);       % vector of actions, used throughout
for s = 1 : ns
    if RL==1
        P_nor{s,1} = (1/No_actions(s,1)) * ones(No_actions(s,1),1);
        P{s,1} = P_nor{s,1}; % P_test{s,1};
        Q{s,1} = P{s,1};
        Y{s,1} = 0;
        P_sto{s,1} = zeros(1,No_actions(s,1));
        R{s,1} = 0;
        R_sto{s,1} = 0;
        Y_sto{s,1} = 0;
    elseif AL==1
        ASP{s,1} = asp_sta;
        ASP_sto{s,1} = zeros(1,1);
        ASP_ave{s,1} = 0;
        ASP_ave_sto{s,1}(1,1) = zeros(1,1);
        % Randomly selection the first action
        P_aggre = zeros(No_actions(s,1),1);
        sum_p = 0;  % sum of probability entries
        for j = 1 : No_actions(s,1)
            P_aggre(j,1) = 1/No_actions(s,1) + sum_p;
            sum_p = P_aggre(j,1);		% update of the sum
        end
        action(s,1) = mixed_strategy(P_aggre',No_actions(s,1));
        % action(s,1) = No_actions(s,1); % last action = connect to everybody
        R{s,1} = 0;
        R_sto{s,1} = 0;
        P{s,1} = zeros(No_actions(s,1),1);  %probability of selecting actions
    elseif AP==1
        % Randomly selection the first action
        P_aggre = zeros(No_actions(s,1),1);
        sum_p = 0;  % sum of probability entries
        for j = 1 : No_actions(s,1)
            P_aggre(j,1) = 1/No_actions(s,1) + sum_p;
            sum_p = P_aggre(j,1);		% update of the sum
        end
        action(s,1) = mixed_strategy(P_aggre',No_actions(s,1));
        % action(s,1) = No_actions(s,1); % last action = connect to everybody
        action_vec{s,1}=vertex(No_actions(s,1),action(s,1));
        for j=1:M
            action_vec_hist{s,1}(:,j)=action_vec{s,1};
        end
        P_nor{s,1} = (1/No_actions(s,1)) * ones(No_actions(s,1),1);
        P{s,1} = P_nor{s,1}; % P_test{s,1};
        R{s,1}=0;
        R_sto{s,1} = 0;
        P_sto{s,1} = zeros(1,No_actions(s,1));
    elseif BB==1
        action(s,1) = 1;  % last action = connect to everybody
        action_base(s,1) = action(s,1);
        R{s,1} = 0;
        R_base{s,1} = R{s,1};
        R_sto{s,1} = R{s,1};
        P_nor{s,1} = (1/No_actions(s,1)) * ones(No_actions(s,1),1);
        P{s,1} = P_nor{s,1}; % P_test{s,1};
    elseif MB==1
        action(s,1) = 1;
        action_base(s,1) = action(s,1);
        mode(s,1) = 0;
        % mode_base(s,1) = mode(s,1); % discontent
        R{s,1} = 0;
        R_base{s,1} = R{s,1};
        R_sto{s,1} = R{s,1};
        % equal probability for the first selected action
        P_nor{s,1} = (1/No_actions(s,1)) * ones(No_actions(s,1),1);
        P{s,1} = P_nor{s,1}; % P_test{s,1};
    elseif TE==1
        action(s,1) = 1;
        action_base(s,1) = action(s,1);
        mode(s,1) = 0;              % discontent
        R{s,1} = 0;
        R_base{s,1} = R{s,1};
        R_sto{s,1} = R{s,1};
        % equal probability for the first selected action
        P_nor{s,1} = (1/No_actions(s,1)) * ones(No_actions(s,1),1);
        P{s,1} = P_nor{s,1}; % P_test{s,1};
    end
    
    R_run_ave{gain,run} = 0;
    R_run_ave_sto{gain,run} = zeros(T_stor,1);
    
    if dist_plot==1
        if loc == 1
            distances{s,1} = zeros(No_neig(s,1)-1,1);
            distances_inf{s,1} = ones(No_neig(s,1)-1,1)*inf;
        elseif loc == 0
            distances{s,1} = zeros(ns,1);
            distances_inf{s,1} = ones(ns,1)*inf;
        end
        %distances_sto{s,1} = zeros(No_neig(s,1)-1,1);
        distances_ave{s,1} = 0;
        distances_ave_sto{s,1} = 0;
        distances_run_ave{s,1} = 0;
        distances_run_ave_sto{s,1} = 0;
        
        inv_distances_run_ave(s,1) = 0;
        inv_distances_run_ave_sto{s,1} = 0;
        total_inv_distances_run_ave(s,1) = 0;
        total_inv_distances_run_ave_sto(1,1) = 0;
    end    
    
    distances_global{s,1} = zeros(ns,1);
    
    production{s,1} = 0;
    deficiency{s,1} = 0;
    No_down_deficients{s,1} = 0;
    down_nodes{s,1}(1,1) = 0;
    No_down_nodes{s,1}(1,1) = 0;
end


%% LOOP STARTS HERE
flash = 1;
t_stor = cell(ns,1);for s=1:ns; t_stor{s,1}=1;end;

for t = 1 : T
    
    if (dist_plot==1)
        sum_distances_run_ave_sto(t_stor{1,1},1) = 0; 
    end
    
    % Here we compute the running average sum of rewards
    if ( t == t_stor_per * (t_stor{1,1}-1) + 1  )
        sum_r = 0;
        for  s = 1:ns
            sum_r = sum_r + R{s,1}/ns;
        end
        R_run_ave{gain,run} = R_run_ave{gain,run} + mu_plot * (sum_r - R_run_ave{gain,run});
        R_run_ave_sto{gain,run}(t_stor{1,1},1) = R_run_ave{gain,run};
        % here we compute the average value
        % (this is essentially the integral over the running average reward
        % received during the simulation time), or the mean running average
        % performance received so far.
        R_run_ave_per_run{gain,1}(run,1) = R_run_ave_per_run{gain,1}(run,1) + R_run_ave_sto{gain,run}(t_stor{1,1},1)/(T_stor);
    end
    total_inv_distances_run_ave = 0;
    
for s = 1 : ns
    
    %------------------
    % STORING VARIABLES
    %------------------
    % For any node
    if t == t_stor_per * (t_stor{s,1}-1) + 1 
        
        if time_plot == 1
            T_stor_vec(t_stor{s,1},1) = (t_stor{s,1}-1) * t_stor_per + 1;
            if RL==1    % for reinforcement learning
                P_sto{s,1}(t_stor{s,1},:) = P{s,1}';
                Y_sto{s,1}(t_stor{s,1},:) = Y{s,1};
                R_sto{s,1}(t_stor{s,1},1) = R{s,1};
            elseif AL==1 % for aspiration learning
                ASP_sto{s,1}(t_stor{s,1},1) = ASP{s,1};
                ASP_ave{s,1} = ASP_ave{s,1} + 1/(t_stor{s,1}+1) * (ASP{s,1} - ASP_ave{s,1});
                ASP_ave_sto{s,1}(t_stor{s,1},1) = ASP_ave{s,1};
                R_sto{s,1}(t_stor{s,1},1) = R{s,1};
            elseif AP==1 % for adaptive play
                P_sto{s,1}(t_stor{s,1},:) = P{s,1}';
                R_sto{s,1}(t_stor{s,1},1) = R{s,1};
            end        
        end
        
        if dist_plot==1 % measuring distances
            T_stor_vec(t_stor{s,1},1) = (t_stor{s,1}-1) * t_stor_per + 1;
            if loc==1
                % if you only consider distances from neighbors
                %distances_sto{s,1}(:,t_stor{s,1}) = distances{s,1};
                distances_ave{s,1} = sum(distances{s,1})/No_neig(s,1);
                 distances_ave_sto{s,1}(1,t_stor{s,1}) = distances_ave{s,1};
                distances_run_ave{s,1}=distances_run_ave{s,1} + 1/(t_stor{s,1}+1) * ... % mu
                    (sum(distances{s,1})/(No_neig(s,1)-1) - distances_run_ave{s,1});
                distances_run_ave_sto{s,1}(t_stor{s,1},1)=distances_run_ave{s,1};
                
                % to compute inverse sum of distances we need a new distances vector
                % where zero distances have been replaced by inf
                for j=1:No_neig(s,1)-1;
                    if distances{s,1}(j,1) > 0
                        distances_inf{s,1}(j,1) = distances{s,1}(j,1);
                    elseif distances{s,1}(j,1) == 0 & j ~= s
                        % if the original distance is zero and j is not node s
                        distances_inf{s,1}(j,1) = inf;
                    elseif distances{s,1}(j,1) == 0 & j==s
                        distances_inf{s,1}(j,1)=0;
                    end
                end
                if dist_plot == 1
                    inv_distances_run_ave(s,1) = inv_distances_run_ave(s,1) +  mu_plot * ... % 1/(t_stor{s,1}+1) * ... %
                        (1/sum(distances_inf{s,1}) - inv_distances_run_ave(s,1));
                    inv_distances_run_ave_sto{s,1}(t_stor{s,1},1) = inv_distances_run_ave(s,1);
                end
                
            elseif loc==0
                % if you consider distances from everybody
                %distances_sto{s,1}(:,t_stor{s,1}) = distances_global{s,1};
                distances_ave{s,1} = sum(distances_global{s,1})/(ns-1);
                distances_ave_sto{s,1}(1,t_stor{s,1}) = distances_ave{s,1};
                distances_run_ave{s,1}=distances_run_ave{s,1} + 1/(t_stor{s,1}+1) * ... % mu
                    (distances_ave{s,1} - distances_run_ave{s,1});
                distances_run_ave_sto{s,1}(t_stor{s,1},1)=distances_run_ave{s,1};
                for j=1:ns;
                    if distances_global{s,1}(j,1) > 0
                        distances_inf{s,1}(j,1) = distances_global{s,1}(j,1);
                    elseif distances_global{s,1}(j,1)==0 & j ~= s
                        % if the original distance is zero
                        distances_inf{s,1}(j,1) = inf;
                    elseif distances_global{s,1}(j,1)==0 & j==s
                        distances_inf{s,1}(j,1) = 0;
                    end
                end
                if (dist_plot == 1)
                    inv_distances_run_ave(s,1) = inv_distances_run_ave(s,1) + mu_plot * ... % 1/(t_stor{s,1}+1) * ... %
                       (1/sum(distances_inf{s,1}) - inv_distances_run_ave(s,1));
                    inv_distances_run_ave_sto{s,1}(t_stor{s,1},1) = inv_distances_run_ave(s,1);
                end
            end
        end
        
        % Computing average total distances
        if (dist_plot == 1)
            total_inv_distances_run_ave = total_inv_distances_run_ave + inv_distances_run_ave(s,1) / ns;
            total_inv_distances_run_ave_sto(t_stor{s,1},1) = total_inv_distances_run_ave;
        end
        t_stor{s,1} = t_stor{s,1} + 1;
        if No_runs==1 && s==1
            disp '---------------'
            disp 'Iteration '
            t
        end
    end
end % of s

%---------------------------------------------------------------
% REINFORCEMENT LEARNING: SELECTING ACTIONS & COMPUTING REWARDS
%---------------------------------------------------------------
if RL==1
    for s = 1:ns
        [P{s,1}, action(s,1) ] = RL_action ( epsilon, P{s,1}, P_nor{s,1}, P_sel{s,1}, Q{s,1}, g(gain,1), No_actions(s,1), PROJ );
    end % of s
    
%------------------------------------------------------------
% ASPIRATION LEARNING: SELECTING ACTIONS & COMPUTING REWARDS
%------------------------------------------------------------
elseif AL==1
    for s=1:ns
        [P{s,1}, action(s,1) ] = AL_action ( ASP{s,1}, R{s,1}, P{s,1}, p_asp, c_asp, No_actions(s,1), action(s,1) );
    end % of s
    
%----------------------------------------
% ADAPTIVE PLAY: SELECTING ACTIONS
%----------------------------------------
elseif AP==1
    % for adaptive play, we need to compute the reward for every available 
    % action, in order to be able to compute a best reply.
    action_prev = action;
    BR = cell(ns,1);
    for s=1:ns  % for each agent
        action(s,1) = AP_action ( action_prev, No_actions, No_neig, Neighbors, con_symaloc, DIAM, ns, s, ...
            rho, decay, kappa_0, kappa_1, distances_global, No_down_nodes, down_nodes, deficiency, const_def, ...
            P, R, BR, Neigh_payoff, p_ada );
    end % of s
     
%----------------------------------------
% BASELINE-BASED DYNAMICS: SELECTING ACTIONS
%----------------------------------------
elseif ( BB == 1 )
    for s = 1 : ns
        action(s,1) = BB_action ( action, action_base, s, epsilon_bb, No_actions, t );
    end
    
%----------------------------------------
% MOOD-BASED DYNAMICS: SELECTING ACTIONS
%----------------------------------------
elseif ( MB == 1 )
    for s = 1 : ns
        action(s,1) = MB_action ( No_actions, action, action_base, s, lambda_mb, c_mb, t );
    end
    
%----------------------------------------
% TRIAL-AND-ERROR DYNAMICS: SELECTING ACTIONS
%----------------------------------------
elseif ( TE == 1 )
    for s = 1 : ns
        action(s,1) = TE_action ( No_actions, action, action_base, mode, lambda_te, s, t );
    end
end

%----------------------
% Computing Rewards
% reward preliminaries
[No_links,L,distances,distances_global,down_nodes,No_down_nodes,LtN] = ...
    reward_preliminaries(ns,DIAM,No_neig,Neighbors,No_actions,action);

% reward
R = reward(ns,DIAM,rho,decay,L,kappa_0,kappa_1,No_neig,Neighbors,con_symaloc,...,
    distances_global,No_down_nodes,down_nodes,deficiency,const_def,No_links,...
    No_actions,action,P,Neigh_payoff);

% noisy observations
for ( s = 1 : ns )
    random_pert_obs = 2 * rand(1) - 1;  % in [-1,1]
    R{s,1} = R{s,1} + random_pert_obs * size_noise;
end

%-----------------------------------------
% REINFORCEMENT LEARNING: UPDATE RECURSION
%-----------------------------------------
if RL==1
for s=1:ns
    action_vec{s,1} = vertex(No_actions(s,1),action(s,1));  % action vector
    if aver == 1 & t>T_aver
        P{s,1} = P{s,1} + (1/(t^q+1)) * R{s,1} * (action_vec{s,1} - P{s,1});
        Q{s,1} = Q{s,1} + (1/(t^q+1)) * (P{s,1} - Q{s,1});
        Y{s,1} = Y{s,1} + (1/(t^q+1)) * (R{s,1} - Y{s,1});
    elseif ((aver == 0 && AR == 0) | (t <= T_aver) )
        P{s,1} = P{s,1} + mu * R{s,1} * (action_vec{s,1} - P{s,1});
        Q{s,1} = Q{s,1} + mu * (P{s,1} - Q{s,1});
        Y{s,1} = Y{s,1} + mu * (R{s,1} - Y{s,1});
    elseif ((aver == 0 && AR == 1) )    % Positive Reinforcement
        P{s,1} = P{s,1} + mu * R{s,1} * (action_vec{s,1} - P{s,1}) * sigmoid_Karandikar(Y{s,1}-R{s,1},p_rl,c_rl);
        Q{s,1} = Q{s,1} + mu * (P{s,1} - Q{s,1});
        Y{s,1} = Y{s,1} + mu * (R{s,1} - Y{s,1});
    end
end % of s

%--------------------------------------
% ASPIRATION LEARNING: UPDATE RECURSION
%--------------------------------------
elseif AL==1

    ASP = AL_update ( asp_pert, ns, ASP, R, s, zeta, asp_min, asp_max, com_pert, mu_asp );

%----------------------------------------
% ADAPTIVE PLAY: UPDATE RECURSION
%----------------------------------------

elseif AP==1
    % the new history of actions
    for s = 1:ns
        action_vec{s,1} = vertex(No_actions(s,1),action(s,1));  % action vector
        action_vec_hist{s,1}(:,2:M) = action_vec_hist{s,1}(:,1:(M-1));
        action_vec_hist{s,1}(:,1) = action_vec{s,1}';
    end
    % given the new history of actions, we compute the new state
    for s = 1:ns
        P{s,1} = sum(action_vec_hist{s,1}')/M;
        P{s,1} = P{s,1}';
    end

%----------------------------------------
% BENCHMARK-BASED: UPDATE BASELINE ACTION
%----------------------------------------
elseif BB==1
    for s = 1 : ns
        action_vec{s,1} = vertex(No_actions(s,1),action(s,1));
        if ( action(s,1) ~= action_base(s,1) )
            % if the action currently selected is different than the base
            % action, then
            if ( R{s,1} > R_base{s,1} )
                action_base(s,1) = action(s,1);
                R_base{s,1} = R{s,1};
            end
            % in any other case, baseline payoff/action remains the same
        else % if it did not experiment
            action_base(s,1) = action(s,1);
            R_base{s,1} = R{s,1};
        end
    end
    
%----------------------------------------
% MOOD-BASED: UPDATE BASELINE ACTION
%----------------------------------------    
elseif MB == 1
    for s = 1 : ns
        [mode(s,1),action_base(s,1),R_base{s,1}] = ...
            MB_update(R{s,1},R_base{s,1},action(s,1),action_base(s,1),mode(s,1),lambda_mb);
    end
   
%----------------------------------------
% TRIAL-AND-ERROR: UPDATE BASELINE ACTION
%----------------------------------------       
elseif TE == 1
    for s = 1 : ns
        % action_vec{s,1} = vertex(No_actions(s,1),action(s,1));
        [mode(s,1),action_base(s,1),R_base{s,1}] = ...
            TE_update(R{s,1},R_base{s,1},action(s,1),action_base(s,1),mode(s,1));
    end
    
end  % of if on Learning Algorithm

%---------------
% PLOTTING LINKS
%---------------
if t == flash_interval * flash & camera == 1
if graph_plot == 1
for s = 1 : ns
    figure(5*((gain-1)*No_runs+run));
    plot(x_c_in{s,1},y_c_in{s,1}),fill(x_c_in{s,1},y_c_in{s,1},'k'),...
        axis([0 1 0 1]),...
        text(x(s,1)+r_in-0.01,y(s,1)+r_in+0.04,num2str(s),'FontSize',10),...
        hold on;
end

for s = 1 : ns
    if No_links(s,1)~=0
        for i = 1 : No_links(s,1)
           arrow([x(LtN{s,1}(i,1),1),y(LtN{s,1}(i,1),1)],[x(s,1),y(s,1)],...
                'tipangle',arrow_tipangle,'length',arrow_length,'Edgecolor','k','Facecolor','k');
            hold on;
        end
    end
end

%------------------
%if camera == 1
    drawnow
    for kk=1:1
        Ff = getframe(gcf,rect);
        avio = addframe(avio,Ff);
    end
    hold off;
    flash = flash + 1;
%end

end
end

end % of t
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%------------
% AVERAGE REWARDS
for s = 1 : ns
    R_ave{gain,1}(run,1) = R_ave{gain,1}(run,1) + R{s,1}/ns;
    D_ave{gain,1}(run,1) = D_ave{gain,1}(run,1) + No_links(s,1)/ns;
end

%--------------------------------------------------------------------------
% Plotting the final actions as a graph
if graph_plot == 1
figure(5*((gain-1)*No_runs+run));
for s = 1 : ns
    figure(5*((gain-1)*No_runs+run));
    plot(x_c_in{s,1},y_c_in{s,1}),fill(x_c_in{s,1},y_c_in{s,1},'k'),...
        axis([0 1 0 1]),...
        text(x(s,1)+r_in+0.01,y(s,1)+r_in+0.04,num2str(s),'FontSize',10),...
        title('Final Graph','FontSize',12),...
        hold on;
end

for s = 1 : ns
    if No_links(s,1)~=0
        for i = 1 : No_links(s,1)
           arrow([x(LtN{s,1}(i,1),1),y(LtN{s,1}(i,1),1)],[x(s,1),y(s,1)],...
                'tipangle',arrow_tipangle,'length',arrow_length,'Edgecolor','k','Facecolor','k');
            hold on;
        end
    end
end

%text(0.1,0.9,['Network value/n = ',num2str(R_ave{gain,1}(run,1),'%6.2f')],...
%    'FontSize',13), hold off;

if camera == 1
    drawnow
    for kk=1:1
        Ff = getframe(gcf,rect);
        avio = addframe(avio,Ff);
    end
    hold off;
    flash = flash + 1;
end

end


%--------------------------------------------------------------------------
% Plotting the responses with time
if time_plot == 1 & ( RL==1 | AP==1 )

%--------------- Plotting the Strategies------------------
if (strategies_plot == 1)
    figure(6*((gain-1)*No_runs+run)+1);
    for s = 1 : ns
        subplot(ns,1,s),
        for j = 1 : No_actions(s,1)
            plot(T_stor_vec(:,1),P_sto{s,1}(1:T_stor,j),char(line_format_col(j))), ...
                hold on;
        end

        ylabel(['Agent ',num2str(s)]),%...
        %axis([0 T -0.1 1.1]),hold on;
    end     % of i
        xlabel('Time step t','FontSize',12)
        legend(char(legend_array(1:No_actions(s,1))),1);
    %     legend(char(legend_array(1:No_actions(s,1))),1,'Location','South',...
    %         'Orientation','horizontal');
        hold off;   
end

%--------------- Plotting the Reward State----------------
% figure(2*((gain-1)*No_runs+run)+5);
% for s = 1 : ns
%     subplot(ns,1,s),
%     for j = 1 : No_actions(s,1)
%         plot(T_stor_vec(:,1),Y_sto{s,1}(1:T_stor,1)), hold on;
%     end
%     xlabel('Time step k'),ylabel(['Agent ',num2str(s)]),...
%     axis([0 T -0.1 1.1]),hold on;
% end     

elseif time_plot==1 & AL==1
    
%--------------- Plotting the Aspirations------------------

% figure(4*((gain-1)*No_runs+run)+2);
% title('Running average aspiration levels vs time')
% for s = 1 : ns
%         subplot(ns,1,s),
%     plot(T_stor_vec(:,1),ASP_ave_sto{s,1}(1:T_stor,1),'k'),
%             axis([0 T 2.5 3]);
%             grid on;hold on;
%     ylabel(['Agent ',num2str(s)]),...
%     %axis([0 T -0.1 1.1]),
%     hold on;
% end     % of i
%     xlabel('Time step k'),
% %     legend(char(legend_array(1:No_actions(s,1))),1,'Location','South',...
% %         'Orientation','horizontal');
%     hold off;   

figure(6*((gain-1)*No_runs+run)+2);
for s = 1 : ns
        subplot(ns,1,s),
    plot(T_stor_vec(:,1),ASP_sto{s,1}(1:T_stor,1),'k'),
            %axis([0 T 1 asp_sta]);
            grid on;hold on;
    ylabel(['Agent ',num2str(s)],'FontSize',12),...
        if s==1; title('Aspiration Level vs Time','FontSize',12); end;
    %axis([0 T -0.1 1.1]),
    hold on;
end     % of i
    xlabel('Time step','FontSize',12),
%     legend(char(legend_array(1:No_actions(s,1))),1,'Location','South',...
%         'Orientation','horizontal');
    hold off;   

end % of plotting if
    
%--------------- Plotting the Distances ------------------
if dist_plot==1
%if ns <=5    
figure(6*((gain-1)*No_runs+run)+3);
for s = 1 : ns
    plot(T_stor_vec(:,1),distances_run_ave_sto{s,1}(1:T_stor,1),...
        char(line_format_col(s))), grid on,
    hold on;
end
legend(char(agent_legend_array(1:ns)),1);
ylabel('Running average of mean distance from neighbors','FontSize',12); 
xlabel('Time step','FontSize',12);
hold off;

%end

% Plotting average distances for all agents
figure(6*((gain-1)*No_runs+run)+4);
for s = 1 : ns
    plot(T_stor_vec(:,1),inv_distances_run_ave_sto{s,1}(1:T_stor,1),char(line_format_col(s))), grid on;
    hold on;
end
legend(char(agent_legend_array(1:ns)),1);
ylabel('Running average of inverse total distance','FontSize',12); 
xlabel('Time step','FontSize',12);
hold off;

% Plotting average distances for all agents
figure(6*((gain-1)*No_runs+run)+5);
% plot(T_stor_vec(:,1),total_inv_distances_run_ave_sto(1:T_stor,1),'k','LineWidth',2), grid on;
plot(T_stor_vec(:,1),R_run_ave_sto{gain,run} ( 1 : T_stor, 1),'k','LineWidth',2), grid on;
% legend(char(agent_legend_array(1:ns)),1);
ylabel('Running average of inverse total distance','FontSize',12); 
xlabel('Time step','FontSize',12);

end % of if

% Note that the inverse distance from the neighbors is a indication of the
% communication energy. The further a node is the smaller number of links
% we use, thus the less energy we use.
% Thus, we can measure the efficiency of the graph by the smaller the
% inverse distance is.
SumPayoffs = sum ( R_run_ave_sto{gain,run} ( 1 : T_stor , 1 ) ) / T_stor ;
SumPayoffs

%------------------------
% Close the camera
if camera == 1
    avio = close(avio);
end

end % of run
end % of gain

%--------------------------------------------------------------------------
if aver_plot == 1 %& RL==1
    % average payoff and degree over all runs
    R_ave_run = zeros(No_gains,No_runs);
    R_run_ave_per_gain = zeros(No_gains,1);
    for gain = 1 : No_gains
        R_ave_run(gain,run) = R_ave{gain,1}(1,1);
        D_ave_run(gain,run) = D_ave{gain,1}(1,1);
        R_run_ave_per_gain(gain,1) = R_run_ave_per_run{gain,1}(1,1);
        for run = 1 : No_runs
    %         R_ave_run(gain,run) = R_ave_run(gain,run) + 0.1 * ( sum(R_ave{gain,1}(1:run,1)) - R_ave_run(gain,run) );
    %             sum(R_ave{gain,1}(1:run,1))/run;
            R_ave_run(gain,run) = sum(R_ave{gain,1}(1:run,1))/run;
            R_run_ave_per_gain(gain,1) = sum(R_run_ave_per_run{gain,1}(1:run,1))/run;
        end
    end

%     figure(1);
%     for gain = 1 : No_gains
%         plot(1:No_runs,R_ave{gain,1},char(line_format_col(gain))),...
%             legend(char(gain_legend_array(1:No_gains)),1),%,'Location','NorthEast'), ...
%             ylabel('Final network value'),xlabel('Run'),hold on;
%     end
% 
%     figure(2);
%     for gain = 1 : No_gains
%         plot(1:No_runs,R_ave_run(gain,1:No_runs),char(line_format_col(gain))),...
%             legend(char(gain_legend_array(1:No_gains)),1),%,'Location','NorthEast'), ...
%             ylabel('Running average final network value per iteration'),xlabel('Run'),hold on;
%     end

    figure(3);
    for gain = 1 : No_gains
        plot(1:No_runs,R_run_ave_per_run{gain,1},char(line_format_col(gain))),...
            legend(char(gain_legend_array(1:No_gains)),1),%,'Location','NorthEast'), ...
            ylabel('Average network value'),xlabel('Run'),hold on;
    end

    figure(4);
    plot(g(1:No_gains),R_run_ave_per_gain(:,1),char(line_format_col(gain))),...
        % legend(char(gain_legend_array(1:No_gains)),1),%,'Location','NorthEast'), ...
        ylabel('Running average payoff per gain'),xlabel('Gain'),hold on;

end

