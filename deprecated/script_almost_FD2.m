if ~exist("conAssign")
    run C:\tomlab\startup.m
else
    disp("TomLab Initiated");
end
addpath(genpath(pwd));
run gen_param.m
%% Initialize estimators
estimator_list = {'FD','AFD','AFD2','HM','EE'};
statistic_list = {'average','bias','var'};
norm_p=[];
norm_p_modified=[];

for estimator = estimator_list
    for statistic = statistic_list
        eval([statistic{1} '_' estimator{1} '=[];']);
    end
end

for nM = [100]
    param.nM = nM;
    [F_struct,state] = DDCMixture.statetransition(param);
    
    F_0 = kron([0,1;0,1],F_struct{2});
    F_1 = kron([1,0;1,0],F_struct{1});

    F   = [F_0;F_1];
    
    F_til = F_1 - F_0;
    f = @(x) obj(x,F_til,F_0);
    Prob = conAssign(f, [], [], [], zeros(64,1), ones(64,1), "Example Problem", zeros(64,1));
    Result = tomRun('snopt', Prob);
    p_star = Result.x_k;
    
    param.P     = F_struct;
    param.state = state;
    param.n_type = 1;
    param.n_action = 2;
%     param.n_state = size(param.P{1},2);
    
%     S1.N     = param.nM;
%     S1.T     = param.nT;
%     S1.beta  = param.beta;
%     S1.P     = F_struct;
%     S1.state = state;
%     S1.n_type = 2;
%     S1.n_action = 2;
%     S1.n_state = size(S1.P{1},2);
    theta_vec = [theta.VP0,theta.VP1,theta.VP2,theta.FC0,theta.FC1,theta.EC0,theta.EC1]';
    ts = tic;
    for i = 1: param.nMC
        [datasim.at,datasim.yt,datasim.zt] = ...
            DDCMixture.simdata(theta_vec,S1,param.nT,param.nM);
        Data{i} = datasim;
    end
    TimeSimulation = toc(ts);
    fprintf('Simulation %d observations of mixture data used %f seconds \n', param.nMC ,TimeSimulation);
    
    norm_p = [norm_p, [Result.f_0]];
    norm_p_modified = [norm_p_modified,[Result.f_k]];
    
    theta_vec0 = zeros(7,1);
    p_default = zeros(64,1);

    % for i = 1:param.nMC
    % Method = {'EE','HM','FD','AFD2'}
    parfor i = 1:param.nMC
        opt = struct();
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
        ResultTable_AFD(i,:) = theta_hat;
        IterTable_AFD(i) = iter;
        TimeTable_AFD(i) = TimeEstimation;    
       

        % The two step AFD with error correctoin
        opt.method = 'AFD2';
        ts = tic;
        [theta_hat,iter] = DDCMixture.SingleEstimation(datasim,S1,theta_vec0,p_star,opt);
        TimeEstimation =  toc(ts);
        ResultTable_AFD2(i,:) = theta_hat;
        IterTable_AFD2(i) = iter;
        TimeTable_AFD2(i) = TimeEstimation;    

    end

    % Put into summary
    for estimator = estimator_list
        eval(['average_' estimator{1} '=[average_' estimator{1}  ' transpose(mean(abs(ResultTable_' estimator{1} ' - transpose(theta_vec))))]; ']);
        eval(['bias_' estimator{1} '=[bias_' estimator{1}  ' transpose(abs(mean(ResultTable_' estimator{1} ' - transpose(theta_vec))))]; ']);
        eval(['var_' estimator{1} '=[var_' estimator{1}  ' transpose(var(ResultTable_' estimator{1} ' - transpose(theta_vec)))]; ']);
    end
end
%% Diary
diarystr = sprintf('diary/AFD_test2_%d_M%d_T%d.txt',param.nGrid,param.nM,param.nT);
delete(diarystr);
diary(diarystr);


% Bias and Variance Table
input.tableRowLabels = {'$\theta_0^{VP}$', '$\theta_1^{VP}$', ...
    '$\theta_2^{VP}$','$\theta_0^{FC}$','$\theta_1^{FC}$','$\theta_0^{EC}$','$\theta_1^{EC}$'};
input.tableColLabels = {'0','0.5','1.0','1.5','2.0'};


for estimator = estimator_list
    eval(['input.data  = bias_' estimator{1} ';']);
    eval(['input.variance= var_' estimator{1} ';']);
    input.tableCaption = ['The bias and variance of ' estimator{1} ' estimator'];
    latexVarianceTable(input);
end

average_time = [mean(TimeTable_HM),mean(TimeTable_EE),mean(TimeTable_FD),mean(TimeTable_AFD),mean(TimeTable_AFD2);
                mean(IterTable_HM),mean(IterTable_EE),mean(IterTable_FD),mean(IterTable_AFD),mean(IterTable_AFD2)];
input.data = average_time;
input.tableColLabels = {'HM','EE','FD','AFD','AFD2'};
input.tableRowLabels = {'Average Time','Average Iter'};
input.tableCaption = 'The averaged time used';
latexTable(input);
diary off;

 
diary off;  