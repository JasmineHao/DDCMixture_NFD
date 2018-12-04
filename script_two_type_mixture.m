runInitEV = 1;
%%
param.nMC = 50;
param.nT = 120;
param.nM = 50;	
param.nGrid = 2;  %number of states for each z_j
param.beta=0.95; % beta used
param.n_action = 2;
param.a_space = [1,0];
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

%parameters that determine the mixture
mix_param.nType = 2;
mix_param.mixProb = [0.5,0.5];

% Run tomLab
if ~exist("conAssign")
    run C:\tomlab\startup.m
else
    disp("TomLab Initiated");
end

diarystr = sprintf('diary/output%d_M%d_T%d.txt',param.nGrid,param.nM,param.nT);
rand_seed = 100;
rand('seed',rand_seed);


%%
%The exogenous transition densities are the same
[P,state] = DDCMixture.statetransition(param);
%Assume the first stage parameters are known
% In this example, S1 and S2 different in the productivity
S1.N     = param.nM;
S1.T     = param.nT;
S1.beta  = param.beta;
S1.P     = P;
S1.state = state;
S1.n_type = 2;
S1.n_action = 2;
S1.n_state = size(P{1},2);
theta_vec = [theta.VP0,theta.VP1,theta.VP2,theta.FC0,theta.FC1,theta.EC0,theta.EC1]';
S2 = S1;
S2.state(:,5) = S1.state(:,5) - 1; 
%--------------------------------------------------------------------------
%Generate Data
%--------------------------------------------------------------------------
%%
ResultTable_NFXP = zeros(param.nMC,7);
TimeTable_NFXP   = zeros(param.nMC,1);
IterTable_NFXP   = zeros(param.nMC,1);
ResultTable_SEQ   = zeros(param.nMC,7);
TimeTable_SEQ     = zeros(param.nMC,1);
IterTable_SEQ   = zeros(param.nMC,1);
ResultTable_FD   = zeros(param.nMC,7);
TimeTable_FD     = zeros(param.nMC,1);
IterTable_FD   = zeros(param.nMC,1);
ResultTable_EE   = zeros(param.nMC,7);
TimeTable_EE     = zeros(param.nMC,1);
IterTable_EE   = zeros(param.nMC,1);

ResultTable_SEQ_EE   = zeros(param.nMC,7);
TimeTable_SEQ_EE     = zeros(param.nMC,1);
IterTable_SEQ_EE   = zeros(param.nMC,1);


ts = tic;
for i = 1: param.nMC
    [datasim.at,datasim.yt,datasim.zt] = ...
        DDCMixture.simdata_mix(S1,S2,theta_vec,param,mix_param);
    Data{i} = datasim;
end
TimeSimulation = toc(ts);
fprintf('Simulation of mixture data used %f seconds \n', TimeSimulation);
%%
%--------------------------------------------------------------------------
%Estimation in NFXP
%--------------------------------------------------------------------------
theta_vec0 = zeros(7,1);
for i = 1:param.nMC
    fprintf('Estimating sample %d out of %d\n', i, param.nMC);
    datasim = Data{i};
    %     ************** NFXP Estimation *****************************
%     ts = tic;
%     opt.method = 'NFXP';
%     [theta_hat,w,iter] = DDCMixture.SequentialEstimation(datasim,mix_param,S1,S2,theta_vec0,opt);
%     TimeEstimation = toc(ts);
%     
%     ResultTable_NFXP(i,:) = theta_hat ;
%     TimeTable_NFXP(i) = TimeEstimation;
%     IterTable_NFXP(i) = iter;
    %     ************** EE Estimation *****************************
%     ts = tic;
%     opt.method = 'EE';
%     [theta_hat,w,iter] = DDCMixture.SequentialEstimation(datasim,mix_param,S1,S2,theta_vec0,opt);
%     ResultTable_EE(i,:) = theta_hat ;
%     TimeEstimation = toc(ts);
%     TimeTable_EE(i) = TimeEstimation;
%     IterTable_EE(i) = iter;
    %     ************** FD Estimation *****************************
    ts = tic;
    opt.method = 'FD';
    [theta_hat,w,iter] = DDCMixture.SequentialEstimation(datasim,mix_param,S1,S2,theta_vec0,opt);
    TimeEstimation =  toc(ts);
    ResultTable_FD(i,:) = theta_hat;
    IterTable_FD(i) = iter;
    TimeTable_FD(i) = TimeEstimation;
    %     ************** EE Estimation *****************************
    ts = tic;
    opt.method = 'SEQ';
    [theta_hat,w,iter] = DDCMixture.SequentialEstimation(datasim,mix_param,S1,S2,theta_vec0,opt);
    TimeEstimation =  toc(ts);
    ResultTable_SEQ(i,:) = theta_hat;
    IterTable_SEQ(i) = iter;
    TimeTable_SEQ(i) = TimeEstimation;
    
    ts = tic;
    opt.method = 'SEQ_EE';
    [theta_hat,w,iter] = DDCMixture.SequentialEstimation(datasim,mix_param,S1,S2,theta_vec0,opt);
    TimeEstimation =  toc(ts);
    ResultTable_SEQ_EE(i,:) = theta_hat;
    IterTable_SEQ_EE(i) = iter;
    TimeTable_SEQ_EE(i) = TimeEstimation;
    
end

%%
diary(diarystr);
fprintf('%d simulations\n',param.nMC);
fprintf('Number of state space : %d\n',S1.n_state * S1.n_action);

average_result = [theta_vec';mean(ResultTable_NFXP);std(ResultTable_NFXP);
    mean(ResultTable_EE);std(ResultTable_EE);
    mean(ResultTable_FD);std(ResultTable_FD);
    mean(ResultTable_SEQ);std(ResultTable_SEQ);
    mean(ResultTable_SEQ_EE);std(ResultTable_SEQ_EE)];
input.data = average_result;
input.dataFormat = {'%.5f'}; 
input.tableColLabels = {'$\theta_0^{VP}$', '$\theta_1^{VP}$', ...
    '$\theta_2^{VP}$','$\theta_0^{FC}$','$\theta_1^{FC}$','$\theta_0^{EC}$','$\theta_1^{EC}$'};
% % Set row labels (use empty string for no label):
input.tableRowLabels = {'DGP','NFXP','','EE','','FD','','SEQ','','SEQ-EE',''}; 
input.tableBorders = 0;
input.tableCaption = 'Simulation results for the market Entry Exit Problem(Two type mixture)';
latexTable(input);

%%

average_result = [mean(abs(ResultTable_NFXP - theta_vec'));
    mean(abs(ResultTable_EE - theta_vec'));
    mean(abs(ResultTable_FD - theta_vec'));
    mean(abs(ResultTable_SEQ - theta_vec'));
    mean(abs(ResultTable_SEQ_EE - theta_vec'))];
input.data = average_result;
input.dataFormat = {'%.5f'}; 
input.tableColLabels = {'$\theta_0^{VP}$', '$\theta_1^{VP}$', ...
    '$\theta_2^{VP}$','$\theta_0^{FC}$','$\theta_1^{FC}$','$\theta_0^{EC}$','$\theta_1^{EC}$'};
% % Set row labels (use empty string for no label):
input.tableRowLabels = {'NFXP','EE','FD','SEQ','SEQ-EE'}; 
input.tableBorders = 0;
input.tableCaption = 'Average Deviation of Entry Exit Problem(Two type mixture)';
latexTable(input);

%%
average_time = [mean(TimeTable_NFXP),mean(TimeTable_EE),mean(TimeTable_FD),mean(TimeTable_SEQ),mean(TimeTable_SEQ_EE)];
average_iter = [mean(IterTable_NFXP),mean(IterTable_EE),mean(IterTable_FD),mean(IterTable_SEQ),mean(IterTable_SEQ_EE)];
input.data = [average_time;average_iter];
input.tableColLabels = {'NFXP','EE','FD','SEQ','SEQ-EE'}; 
input.tableRowLabels = {'time','iter'};

%
latexTable(input);


diary off;
%