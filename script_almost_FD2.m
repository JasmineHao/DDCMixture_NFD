param.nMC = 500;
param.nT = 120;
param.nM = 50;	
param.nGrid = 2;  %number of states for each z_j
param.beta=0.95; % beta used
param.n_action = 2;
param.a_space = [1,0];
% Exogeneous state variable that determine the exogenous transition
param.min_z = 0;
param.max_z = 1;

param.gamma.gamma_a = 5;
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

theta.VP0 = 1.5;
theta.VP1 = 1.0;
theta.VP2 = -1.0;
theta.FC0 = 0.5;
theta.FC1 = 1.0;
theta.EC0 = 1.0;
theta.EC1 = 1.0;
theta.pnames = {'VP0','VP1','VP2','FC0','FC1','EC0','EC1'};

if ~exist("conAssign")
    run C:\tomlab\startup.m
else
    disp("TomLab Initiated");
end

%%
average_FD          = [];
average_FD_modified = [];
bias_FD             = [];
bias_FD_modified    = [];
var_FD              = [];
var_FD_modified     = [];
norm_p              = [];
norm_p_modified     = [];


average_HM          = [];
average_EE          = [];
bias_HM             = [];
bias_EE             = [];
var_HM              = [];
var_EE              = [];


for gamma_a = 0 : 0.5 : 2.0
% for gamma_a = 0.5 : 0.5 
    
    param.gamma.gamma_a = gamma_a;
    [F_struct,state] = DDCMixture.statetransition(param);
    
    F_0 = kron([0,1;0,1],F_struct{2});
    F_1 = kron([1,0;1,0],F_struct{1});

    F   = [F_0;F_1];
    
    F_til = F_1 - F_0;
    f = @(x) obj(x,F_til,F_0);
    Prob = conAssign(f, [], [], [], zeros(64,1), ones(64,1), "Example Problem", zeros(64,1));
    Result = tomRun('snopt', Prob);
    p_star = Result.x_k;
    
    
    S1.N     = param.nM;
    S1.T     = param.nT;
    S1.beta  = param.beta;
    S1.P     = F_struct;
    S1.state = state;
    S1.n_type = 2;
    S1.n_action = 2;
    S1.n_state = size(S1.P{1},2);
    theta_vec = [theta.VP0,theta.VP1,theta.VP2,theta.FC0,theta.FC1,theta.EC0,theta.EC1]';
    ts = tic;
    for i = 1: param.nMC
        [datasim.at,datasim.yt,datasim.zt] = ...
            DDCMixture.simdata(theta_vec,S1,param.nT,param.nM);
        Data{i} = datasim;
    end
    TimeSimulation = toc(ts);
    fprintf('Simulation %d observations of mixture data used %f seconds \n', param.nMC ,TimeSimulation);
    %     disp(Result.x_k);
    %     disp(Result.f_k);
    %     disp(Result.f_0);
    
    ResultTable_FD   = zeros(param.nMC,7);
    TimeTable_FD     = zeros(param.nMC,1);
    IterTable_FD     = zeros(param.nMC,1);

    ResultTable_FD_modified   = zeros(param.nMC,7);
    TimeTable_FD_modified     = zeros(param.nMC,1);
    IterTable_FD_modified     = zeros(param.nMC,1);


    ResultTable_EE   = zeros(param.nMC,7);
    TimeTable_EE     = zeros(param.nMC,1);
    IterTable_EE     = zeros(param.nMC,1);

    ResultTable_HM_modified   = zeros(param.nMC,7);
    TimeTable_HM_modified     = zeros(param.nMC,1);
    IterTable_HM_modified     = zeros(param.nMC,1);
    
    norm_p = [norm_p, [Result.f_0]];
    norm_p_modified = [norm_p_modified,[Result.f_k]];
    
    theta_vec0 = zeros(7,1);
    p_default = zeros(64,1);

    % for i = 1:param.nMC
    for i = 1:param.nMC
        fprintf('Estimating sample %d out of %d\n', i, param.nMC);
        datasim = Data{i};
        
        ts = tic;
        opt.method = 'EE';
        [theta_hat,iter] = DDCMixture.SingleEstimation(datasim,S1,theta_vec0,p_star,opt);
        TimeEstimation =  toc(ts);
        ResultTable_EE(i,:) = theta_hat;
        IterTable_EE(i) = iter;
        TimeTable_EE(i) = TimeEstimation;    

        ts = tic;
        opt.method = 'HM';
        [theta_hat,iter] = DDCMixture.SingleEstimation(datasim,S1,theta_vec0,p_star,opt);
        TimeEstimation =  toc(ts);
        ResultTable_HM(i,:) = theta_hat;
        IterTable_HM(i) = iter;
        TimeTable_HM(i) = TimeEstimation;    

        
        ts = tic;
        opt.method = 'FD';
        [theta_hat,iter] = DDCMixture.SingleEstimation(datasim,S1,theta_vec0,p_default,opt);
        TimeEstimation =  toc(ts);
        ResultTable_FD(i,:) = theta_hat;
        IterTable_FD(i) = iter;
        TimeTable_FD(i) = TimeEstimation;    


        ts = tic;
        [theta_hat,iter] = DDCMixture.SingleEstimation(datasim,S1,theta_vec0,p_star,opt);
        TimeEstimation =  toc(ts);
        ResultTable_FD_modified(i,:) = theta_hat;
        IterTable_FD_modified(i) = iter;
        TimeTable_FD_modified(i) = TimeEstimation;    
       
        
    end

    average_FD = [average_FD  mean(abs(ResultTable_FD - theta_vec'))'];
    average_FD_modified = [average_FD_modified  mean(abs(ResultTable_FD_modified - theta_vec'))'];
    average_EE = [average_EE  mean(abs(ResultTable_EE - theta_vec'))'];
    average_HM = [average_HM  mean(abs(ResultTable_HM - theta_vec'))'];
    bias_FD = [bias_FD abs(mean(ResultTable_FD - theta_vec'))'];
    bias_FD_modified = [bias_FD_modified  abs(mean(ResultTable_FD_modified - theta_vec'))'];
    bias_EE = [bias_EE abs(mean(ResultTable_EE - theta_vec'))'];
    bias_HM = [bias_HM abs(mean(ResultTable_HM - theta_vec'))'];
    var_FD = [var_FD  var(ResultTable_FD )'];
    var_FD_modified = [var_FD_modified  var(ResultTable_FD_modified)'];
    var_EE = [var_EE  var(ResultTable_EE )'];
    var_HM = [var_HM  var(ResultTable_HM )'];
end
%%
diarystr = sprintf('diary/AFD%d_M%d_T%d.txt',param.nGrid,param.nM,param.nT);
diary(diarystr);

input.data = [0 : 0.5 : 2.0; norm_p;norm_p_modified];
input.tableRowLabels = {'$\gamma_a$', 'norm before modified', 'norm after modified'}; 
input.tableColLabels = { '1','2','3','4','5' };
input.tableCaption = 'The norm of the differences in transition densities';
latexTable(input);


input.tableRowLabels = {'$\theta_0^{VP}$', '$\theta_1^{VP}$', ...
    '$\theta_2^{VP}$','$\theta_0^{FC}$','$\theta_1^{FC}$','$\theta_0^{EC}$','$\theta_1^{EC}$'};
input.tableColLabels =  num2cell(0 : 0.5 : 2.0) ;


input.data = average_FD;
input.tableCaption = 'The average deviation of finite dependence estimator';
latexTable(input);

input.data = average_FD_modified;
input.tableCaption = 'The average deviation of almost finite dependence estimator';
latexTable(input);

input.data = average_EE;
input.tableCaption = 'The average deviation of EE estimator';
latexTable(input);

input.data = average_HM;
input.tableCaption = 'The average deviation of HM dependence estimator';
latexTable(input);



input.data    = bias_FD;
input.variance= var_FD;
input.tableCaption = 'The bias and variance of finite dependence estimator';
latexVarianceTable(input);


input.data    = bias_FD_modified;
input.variance= var_FD_modified;
input.tableCaption = 'The bias and variance of almost finite dependence estimator';
latexVarianceTable(input);

input.data    = bias_EE;
input.variance= var_EE;
input.tableCaption = 'The bias and variance of EE estimator';
latexVarianceTable(input);

input.data    = bias_HM;
input.variance= var_HM;
input.tableCaption = 'The bias and variance of HM estimator';
latexVarianceTable(input);


diary off;

% Prob = conAssign(f, g, H, HessPattern, x_L, x_U, Name, x_0, ...
%                             pSepFunc, fLowBnd, ...
%                             A, b_L, b_U, c, dc, d2c, ConsPattern, c_L, c_U, ...
%                             x_min, x_max, f_opt, x_opt);
 
%% Use tomlab to solve

