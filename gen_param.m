runInitEV = 1;
param.nMC = 5;
param.nT = 120;
param.nM = 50;	
param.nGrid = 2;  %number of states for each z_j
param.beta=0.95; % beta used
param.n_action = 2;
param.a_space = [1,0];
param.n_state = param.nGrid^5;
% Exogeneous state variable that determine the exogenous transition
param.min_z = 0;
param.max_z = 1;
param.gamma.gamma_a = 0.5;
param.min_omega = -1;
param.max_omega = 1;
param.gamma.z_0 = 0; %Parameters of z_j transition
param.gamma.z_1 = 0.6; 
param.gamma.omega_0 = 0.2; % Parameters of omega transition(productivity)
param.gamma.omega_1 = 0.9;
param.gamma.sigma_z = 1;
param.gamma.sigma_omega = 1;

% Two thetas
theta.nparam = 7;
theta.VP0 = 0.5;
theta.VP1 = 1.0;
theta.VP2 = -1.0;
theta.FC0 = 0.5;
theta.FC1 = 1.0;
theta.EC0 = 1.0;
theta.EC1 = 1.0;
theta.pnames = {'VP0','VP1','VP2','FC0','FC1','EC0','EC1'};

%system and MC parameters
param.MC=1;					  %number of MC iterations
param.n_type=1;
param.N=param.nM;
param.T=param.nT;
