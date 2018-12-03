

%parameters that determine the mixture
mix_param.nType = 2;
mix_param.mixProb = [0.5,0.5];

% Run tomLab
if ~exist("conAssign")
    run C:\tomlab\startup.m
else
    disp("TomLab Initiated");
end

diarystr = sprintf('diary/output_M%d_T%d.txt',param.nM,param.nT);
rand_seed = 100;
rand('seed',rand_seed);
%%
%The exogenous transition densities are the same
[P,state] = DDCMixture.statetransition(param);
%Assume the first stage parameters are known
% In this example, S and S2 different in the productivity
S.N     = param.nM;
S.T     = param.nT;
S.beta  = param.beta;
S.P     = P;
S.state = state;
S.n_type = 2;
S.n_action = 2;
S.n_state = size(P{1},2);
theta_vec = [theta.VP0,theta.VP1,theta.VP2,theta.FC0,theta.FC1,theta.EC0,theta.EC1]';

%%
%Generate Data

ts = tic;
for i = 1: param.nMC
    [datasim.at,datasim.yt,datasim.zt] = ...
        DDCMixture.simdata(theta_vec,S,param.nT,param.nM);
    Data{i} = datasim;
end
TimeSimulation = toc(ts);
fprintf('Simulation of mixture data used %f seconds \n', TimeSimulation);
