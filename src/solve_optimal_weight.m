function w=solve_optimal_weight(param,F_struct)
    x_size=param.n_action*param.n_state;
    d_size=param.n_action-1;
    row_index=[];
    for k = 1:x_size
        for m = 1:d_size
            row_index = [row_index,(k-1)*(x_size*d_size) + (m-1)*x_size+k];
        end
    end

    %
%     [F_struct,state] = DDCMixture.statetransition(param); %Generate
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
    w = inv(A'*A) * A'* (-b);
end