clear all;
% if ~exist("conAssign")
%     run C:\tomlab\startup.m
% else
%     disp("TomLab Initiated");
% end
addpath(genpath(pwd));
run gen_param.m
max_iter=1;
param.nMC=100;

%% Initialize estimators
% This experiment mainly shows that as the 
estimator_list = {'FD','FD2','AFD','AFD2','HM','EE','HM_true','EE_true'};
statistic_list = {'average','bias','var','time','iter'};
gamma_a_list = [0,1,2,5];
norm_p=[];
norm_p_modified=[];

for estimator = estimator_list
    for statistic = statistic_list
        eval([statistic{1} '_' estimator{1} '=[];']);
    end
end

for gamma_a = gamma_a_list
    
    param.gamma.gamma_a = gamma_a;
    [F_struct,state] = DDCMixture.statetransition(param); %Generate
    
    F_0 = kron([0,1;0,1],F_struct{2});
    F_1 = kron([1,0;1,0],F_struct{1});

    F   = [F_0;F_1];
    
    F_til = F_1 - F_0;
    f = @(x) obj(x,F_til,F_0);
    [p_star,f_val] = fmincon(f,zeros(64,1),[],[],[],[],1e-6*ones(64,1),(1-1e-6*ones(64,1)))
    
    
    %Prob = conAssign(f, [], [], [], zeros(64,1), ones(64,1), "Example Problem", zeros(64,1));
    %Result = tomRun('snopt', Prob);
    p_star = Result.x_k;
    
    param.P     = F_struct;
    param.state = state;
    param.n_type = 1;
    param.n_action = 2;
    
    theta_vec = [theta.VP0,theta.VP1,theta.VP2,theta.FC0,theta.FC1,theta.EC0,theta.EC1]';
    ts = tic;
    for i = 1: param.nMC
        [datasim.at,datasim.yt,datasim.zt] = ...
            DDCMixture.simdata(theta_vec,param,param.nT,param.nM);
        Data{i} = datasim;
    end
    
    ev=zeros(param.n_state,param.n_action);
    pi = DDCMixture.dpidth(param) * theta_vec;
    [p1,ev] = DDCMixture.solveNFXP(ev,pi,param); 
    param.p_1=p1;
    
    TimeSimulation = toc(ts);
    fprintf('Simulation %d observations of mixture data used %f seconds \n', param.nMC ,TimeSimulation);
    
    norm_p = [norm_p, [Result.f_0]];
    norm_p_modified = [norm_p_modified,[Result.f_k]];
    
    theta_vec0 = zeros(7,1);
    p_default = zeros(64,1);

 
    parfor i = 1:param.nMC
        opt = struct();
        
        fprintf('Estimating sample %d out of %d\n', i, param.nMC);
        datasim = Data{i};
        
        opt.true_ccp=0;
        ts = tic;
        opt.method = 'EE'; opt.max_iter=max_iter; %The sequential version
        [theta_hat,iter] = DDCMixture.SingleEstimation(datasim,param,theta_vec0,p_star,opt);
        TimeEstimation =  toc(ts);
        ResultTable_EE(i,:) = theta_hat;
        IterTable_EE(i) = iter;
        TimeTable_EE(i) = TimeEstimation;    

        ts = tic;
        opt.method = 'HM';opt.max_iter=max_iter;
        [theta_hat,iter] = DDCMixture.SingleEstimation(datasim,param,theta_vec0,p_star,opt);
        TimeEstimation =  toc(ts);
        ResultTable_HM(i,:) = theta_hat;
        IterTable_HM(i) = iter;
        TimeTable_HM(i) = TimeEstimation;    

        
        ts = tic;
        opt.method = 'FD';opt.max_iter=max_iter;
        [theta_hat,iter] = DDCMixture.SingleEstimation(datasim,param,theta_vec0,p_default,opt);
        TimeEstimation =  toc(ts);
        ResultTable_FD(i,:) = theta_hat;
        IterTable_FD(i) = iter;
        TimeTable_FD(i) = TimeEstimation;    


        ts = tic;
        [theta_hat,iter] = DDCMixture.SingleEstimation(datasim,param,theta_vec0,p_star,opt);
        TimeEstimation =  toc(ts);
        ResultTable_AFD(i,:) = theta_hat;
        IterTable_AFD(i) = iter;
        TimeTable_AFD(i) = TimeEstimation;    
       

%         The two step AFD with error correctoin
        
        opt.method = 'FD2';
        ts = tic;opt.max_iter=max_iter;
        [theta_hat,iter] = DDCMixture.SingleEstimation(datasim,param,theta_vec0,p_star,opt);
        TimeEstimation =  toc(ts);
        ResultTable_AFD2(i,:) = theta_hat;
        IterTable_AFD2(i) = iter;
        TimeTable_AFD2(i) = TimeEstimation;    

        opt.method = 'FD2';
        ts = tic;opt.max_iter=max_iter;
        [theta_hat,iter] = DDCMixture.SingleEstimation(datasim,param,theta_vec0,p_default,opt);
        TimeEstimation =  toc(ts);
        ResultTable_FD2(i,:) = theta_hat;
        IterTable_FD2(i) = iter;
        TimeTable_FD2(i) = TimeEstimation;    

        opt.true_ccp=1;
        ts = tic;
        opt.method = 'EE'; opt.max_iter=max_iter; %The sequential version
        [theta_hat,iter] = DDCMixture.SingleEstimation(datasim,param,theta_vec0,p_star,opt);
        TimeEstimation =  toc(ts);
        ResultTable_EE_true(i,:) = theta_hat;
        IterTable_EE_true(i) = iter;
        TimeTable_EE_true(i) = TimeEstimation;    

        ts = tic;
        opt.method = 'HM';opt.max_iter=max_iter;
        [theta_hat,iter] = DDCMixture.SingleEstimation(datasim,param,theta_vec0,p_star,opt);
        TimeEstimation =  toc(ts);
        ResultTable_HM_true(i,:) = theta_hat;
        IterTable_HM_true(i) = iter;
        TimeTable_HM_true(i) = TimeEstimation;    
    end
    % Put into summary
    for estimator = estimator_list
        eval(['average_' estimator{1} '=[average_' estimator{1}  ' transpose(mean(abs(ResultTable_' estimator{1} ' - transpose(theta_vec))))]; ']);
        eval(['bias_' estimator{1} '=[bias_' estimator{1}  ' transpose(abs(mean(ResultTable_' estimator{1} ' - transpose(theta_vec))))]; ']);
        eval(['var_' estimator{1} '=[var_' estimator{1}  ' transpose(var(ResultTable_' estimator{1} ' - transpose(theta_vec)))]; ']);
        eval(['time_' estimator{1} '=[time_' estimator{1}  ' transpose(TimeTable_' estimator{1} ')]; ']);
        eval(['iter_' estimator{1} '=[iter_' estimator{1}  ' transpose(IterTable_' estimator{1} ')]; ']);
    end
end
%% Diary Session
diarystr = sprintf('diary/Table_2step_gammaa_%d_M%d_T%d.txt',param.nGrid,param.nM,param.nT);
delete(diarystr);
diary(diarystr);
disp(['This experiment uses 2 step estimator with different values of'...
    '$\gamma_a$. The size of the sample is N=100, T=120 with 500 Monte Carlo'...
    'Simulations.']);

input.data = [norm_p;norm_p_modified];
input.tableRowLabels = {'norm before modified', 'norm after modified'}; 
input.tableColLabels = num2cell(gamma_a_list);
input.tableCaption = 'The norm of the differences in transition densities';
latexTable(input);


% Bias and Variance Table
input.tableRowLabels = param_cell;
input.tableColLabels = num2cell(gamma_a_list);


for estimator = estimator_list
    eval(['input.data  = bias_' estimator{1} ';']);
    eval(['input.variance= var_' estimator{1} ';']);
    input.tableCaption = ['The bias and variance of ' estimator{1} ' estimator'];
    input.tableLabel=['2step' estimator{1}];
    latexVarianceTable(input);
end


command_str_time = '[';command_str_iter = '[';
for estimator = estimator_list
    command_str_time = [command_str_time '; mean(time_',estimator{1},') '];
    command_str_iter = [command_str_iter '; mean(iter_',estimator{1},') '];
end
command_str_time = [command_str_time ']'];command_str_iter = [command_str_iter ']'];

input.data = eval(command_str_time);
input.variance = eval(command_str_iter);
input.tableColLabels = num2cell(gamma_a_list);
input.tableRowLabels = estimator_list;
input.tableCaption = 'The averaged time used';
latexVarianceTable(input);
diary off;
