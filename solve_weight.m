clear all;
if ~exist("conAssign")
    run C:\tomlab\startup.m
else
    disp("TomLab Initiated");
end
addpath(genpath(pwd));
run gen_param.m
max_iter=1;
%%
gamma_a=1;
param.gamma.gamma_a = gamma_a;
%%
param.gamma.gamma_a =1;
param.nGrid=2;
param.n_state=param.nGrid^5;
[F_struct,state] = DDCMixture.statetransition(param); %Generate
x_size=param.n_action*param.n_state;
d_size=param.n_action-1;
row_index=[];
for k = 1:x_size
    for m = 1:d_size
        row_index = [row_index,(k-1)*(x_size*d_size) + (m-1)*x_size+k];
    end
end
%
[F_struct,state] = DDCMixture.statetransition(param); %Generate
F_0 = kron([0,1;0,1],F_struct{2});
F_1 = kron([1,0;1,0],F_struct{1});
F   = [F_0;F_1];    
F_til = F_1 - F_0;
    
% A = kron(F_til',F_til);
% A = A(:,row_index');
%b = b(row_index');
%
a = F_til';
b = F_til;
nrow_b = size(b,2);
A = [];
for each_ind = row_index
    ind_b = rem(each_ind, nrow_b);
    if ind_b == 0
        ind_b = nrow_b;
        ind_a = floor(each_ind / nrow_b) ;
    else
        ind_a = floor(each_ind / nrow_b) + 1;
    end
    
    A = [A,kron(a(:,ind_a),b(:,ind_b) )];
end
b = -vec(F_til*F_0);
w1 = inv(A'*A) * A'* b;

f = @(x) obj(x,F_til,F_0);
Prob = conAssign(f, [], [], [], zeros(64,1), ones(64,1), "Example Problem", zeros(64,1));
Result = tomRun('snopt', Prob);
w2 = (Result.x_k);
disp(f(w1));
disp(f(w2));
%%
% may suffer memory error(must suffer), need to find more efficient solving
% method.
param.nGrid=3;
param.n_state=param.nGrid^5;
[F_struct,state] = DDCMixture.statetransition(param); %Generate
x_size=param.n_action*param.n_state;
d_size=param.n_action-1;
row_index=[];
for k = 1:x_size
    for m = 1:d_size
        row_index = [row_index,(k-1)*(x_size*d_size) + (m-1)*x_size+k];
    end
end

%
[F_struct,state] = DDCMixture.statetransition(param); %Generate
F_0 = kron([0,1;0,1],F_struct{2});
F_1 = kron([1,0;1,0],F_struct{1});
F   = [F_0;F_1];    
F_til = F_1 - F_0;
    
a = F_til';
b = F_til;
nrow_b = size(b,2);
A = [];

for each_ind = row_index
    ind_b = rem(each_ind, nrow_b);
    if ind_b == 0
        ind_b = nrow_b;
        ind_a = floor(each_ind / nrow_b) ;
    else
        ind_a = floor(each_ind / nrow_b) + 1;
    end
    
    A = [A,kron(a(:,ind_a),b(:,ind_b) )];
end
b = vec(F_til*F_0);
w1 = inv(A'*A) * A'* (-b);
%%
f1 = @(x) norm(A * x + b);
Prob = conAssign(f1, [], [], [], zeros(x_size,1), ones(x_size,1), "Example Problem", zeros(x_size,1));
Result = tomRun('KNITRO', Prob);
w1 = (Result.x_k);
 %%
f = @(x) obj(x,F_til,F_0);
Prob = conAssign(f, [], [], [], zeros(x_size,1), ones(x_size,1), "Example Problem", zeros(x_size,1));
Result = tomRun('snopt', Prob);
w2 = (Result.x_k);