clear all;
addpath(genpath(pwd));
run gen_param.m
max_iter=1000;
param.nMC=50;

%% Initialize estimators
% This experiment mainly shows that as the
estimator_list = {'HM','EE'};
statistic_list = {'average','bias','var','time','iter'};
gamma_a_list = [0,1,2,3,4,5];
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
    
    param.P     = F_struct;
    param.state = state;
    param.n_type = 1;
    param.n_action = 2;

    theta_vec = [theta.VP0,theta.VP1,theta.VP2,theta.FC0,theta.FC1,theta.EC0,theta.EC1]';
    ts=tic;
    for i = 1: param.nMC
        [datasim.at,datasim.yt,datasim.zt] = ...
            DDCMixture.simdata(theta_vec,param,param.nT,param.nM);
        Data{i} = datasim;
    end
    ev=zeros(param.n_state,param.n_action);
    pi = DDCMixture.dpidth(param) * theta_vec;
    [p1,ev] = DDCMixture.solveNFXP(ev,pi,param); 
    param.p1=p1;

    TimeSimulation = toc(ts);
    fprintf('Simulation %d observations of mixture data used %f seconds \n', param.nMC ,TimeSimulation);

    theta_vec0 = zeros(7,1);
    p_default = zeros(param.n_state*param.n_action,1);

%%
    parfor i = 1:param.nMC
        
        opt = struct();
        opt.true_ccp=0;
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


    % Put into summary
        for estimator = estimator_list
            eval(['average_' estimator{1} '=[average_' estimator{1}  ' transpose(mean(abs(ResultTable_' estimator{1} ' - transpose(theta_vec))))]; ']);
            eval(['bias_' estimator{1} '=[bias_' estimator{1}  ' transpose(abs(mean(ResultTable_' estimator{1} ' - transpose(theta_vec))))]; ']);
            eval(['var_' estimator{1} '=[var_' estimator{1}  ' transpose(var(ResultTable_' estimator{1} ' - transpose(theta_vec)))]; ']);
            eval(['time_' estimator{1} '=[time_' estimator{1}  ' transpose(TimeTable_' estimator{1} ')]; ']);
            eval(['iter_' estimator{1} '=[iter_' estimator{1}  ' transpose(IterTable_' estimator{1} ')]; ']);
        end
    end
end