if ~exist("conAssign")
    run C:\tomlab\startup.m
else
    disp("TomLab Initiated");
end
addpath(genpath(pwd));
run gen_param.m
max_iter=1;
%% Initialize estimators
estimator_list = {'FD','FD2','AFD','AFD2','HM','EE'};
statistic_list = {'average','bias','var'};
norm_p=[];
norm_p_modified=[];

for estimator = estimator_list
    for statistic = statistic_list
        eval([statistic{1} '_' estimator{1} '=[];']);
    end
end

for gamma_a = 0 : 0.5 : 2.0
% for gamma_a = 2.0 : 0.5 : 3.0
     param.gamma.gamma_a = gamma_a;
    [F_struct,state] = DDCMixture.statetransition(param); %Generate
    
    F_0 = kron([0,1;0,1],F_struct{2});
    F_1 = kron([1,0;1,0],F_struct{1});

    F   = [F_0;F_1];
    
    F_til = F_1 - F_0;
    f = @(x) obj(x,F_til,F_0);
    Prob = conAssign(f, [], [], [], zeros(64,1), ones(64,1), "Example Problem", zeros(64,1));
    Result = tomRun('snopt', Prob);
    p_star = Result.x_k;
    
    
%     S1.N     = param.nM;
%     S1.T     = param.nT;
%     S1.beta  = param.beta;
    param.P     = F_struct;
    param.state = state;
    param.n_type = 1;
    param.n_action = 2;
%     param.n_state = size(S1.P{1},2);
    
    theta_vec = [theta.VP0,theta.VP1,theta.VP2,theta.FC0,theta.FC1,theta.EC0,theta.EC1]';
    ts = tic;
    for i = 1: param.nMC
        [datasim.at,datasim.yt,datasim.zt] = ...
            DDCMixture.simdata(theta_vec,param,param.nT,param.nM);
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
        
        opt.method = 'AFD2';
        ts = tic;opt.max_iter=max_iter;
        [theta_hat,iter] = DDCMixture.SingleEstimation(datasim,param,theta_vec0,p_star,opt);
        TimeEstimation =  toc(ts);
        ResultTable_AFD2(i,:) = theta_hat;
        IterTable_AFD2(i) = iter;
        TimeTable_AFD2(i) = TimeEstimation;    

        opt.method = 'AFD2';
        ts = tic;opt.max_iter=max_iter;
        [theta_hat,iter] = DDCMixture.SingleEstimation(datasim,param,theta_vec0,p_default,opt);
        TimeEstimation =  toc(ts);
        ResultTable_FD2(i,:) = theta_hat;
        IterTable_FD2(i) = iter;
        TimeTable_FD2(i) = TimeEstimation;    

    end
    % Put into summary
    for estimator = estimator_list
        eval(['average_' estimator{1} '=[average_' estimator{1}  ' transpose(mean(abs(ResultTable_' estimator{1} ' - transpose(theta_vec))))]; ']);
        eval(['bias_' estimator{1} '=[bias_' estimator{1}  ' transpose(abs(mean(ResultTable_' estimator{1} ' - transpose(theta_vec))))]; ']);
        eval(['var_' estimator{1} '=[var_' estimator{1}  ' transpose(var(ResultTable_' estimator{1} ' - transpose(theta_vec)))]; ']);
    end
end
%% Diary Session
diarystr = sprintf('diary/Table_2step_gammaa_%d_M%d_T%d.txt',param.nGrid,param.nM,param.nT);
delete(diarystr);
diary(diarystr);

input.data = [0 : 0.5 : 2.0; norm_p;norm_p_modified];
input.tableRowLabels = {'$\gamma_a$', 'norm before modified', 'norm after modified'}; 
input.tableColLabels = { '','','','','' };
input.tableCaption = 'The norm of the differences in transition densities';
latexTable(input);


% Bias and Variance Table
input.tableRowLabels = param_cell;
input.tableColLabels = {'0','0.5','1.0','1.5','2.0'};


for estimator = estimator_list
    eval(['input.data  = bias_' estimator{1} ';']);
    eval(['input.variance= var_' estimator{1} ';']);
    input.tableCaption = ['The bias and variance of ' estimator{1} ' estimator'];
    input.tableLabel=['2step' estimator{1}];
    latexVarianceTable(input);
end


command_str_time = '[';command_str_iter = '[';
for estimator = estimator_list
    command_str_time = [command_str_time 'mean(TimeTable_',estimator{1},') '];
    command_str_iter = [command_str_iter 'mean(IterTable_',estimator{1},') '];
end
command_str_time = [command_str_time ']'];command_str_iter = [command_str_iter ']'];


average_time = [eval(command_str_time);eval(command_str_iter)];
input.data = average_time;
input.tableColLabels = estimator_list;
input.tableRowLabels = {'Average Time','Average Iter'};
input.tableCaption = 'The averaged time used';
latexTable(input);
diary off;
