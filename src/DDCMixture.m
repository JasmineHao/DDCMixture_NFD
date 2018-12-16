classdef DDCMixture
    methods (Static)
        function [P,state] = statetransition(param)
            %--------------------------------------------------------------
            % NFXP.statetransition
            % Input: param
            % Output:
            %         P:        n_state x n_state vector
            %                   state transition matrix
            %
            %         state:    n_state x dim(state)
            %                   define all possible states
            %--------------------------------------------------------------
            min_z = param.min_z;
            max_z = param.max_z;
            min_omega = param.min_omega;
            max_omega = param.max_omega;
            gam_z0 = param.gamma.z_0; %Parameters of z_j transition
            gam_z1 = param.gamma.z_1;
            gam_omega0 = param.gamma.omega_0; % Parameters of omega transition(productivity)
            gam_omega1 = param.gamma.omega_1;
            sigma_z = param.gamma.sigma_z;
            sigma_omega = param.gamma.sigma_omega;
            gamma_a = param.gamma.gamma_a;
            nGrid = param.nGrid;
            
            w_z = (max_z - min_z)/(nGrid - 1);
            w_omega = (max_omega - min_omega)/(nGrid - 1);
            z_grid = min_z : w_z : max_z;
            o_grid = min_omega : w_omega : max_omega;
            %--------------------------------------------------------------
            % f_z: transition density for z_j
            f_z = zeros(nGrid,nGrid);
            for i = 1 : nGrid
                for j = 1: nGrid
                    if j == 1
                        f_z(i,j) = normcdf((z_grid(j) + w_z / 2 - gam_z0 - ...
                            gam_z1 * z_grid(i))/sqrt(sigma_z));
                    elseif j == nGrid
                        f_z(i,j) = 1 - normcdf((z_grid(j) - w_z / 2 - gam_z0...
                            - gam_z1 * z_grid(i))/sqrt(sigma_z));
                    else
                        f_z(i,j) = normcdf((z_grid(j) + w_z / 2 - gam_z0 - ...
                            gam_z1 * z_grid(i))/sqrt(sigma_z)) - ...
                            normcdf((z_grid(j) - w_z / 2 - gam_z0 - ...
                            gam_z1 * z_grid(i) )/sqrt(sigma_z));
                    end
                end
            end
            %--------------------------------------------------------------
            % f_o: transition density for omega_j
            F_z = kron(f_z,kron(f_z,kron(f_z,f_z)));
            f_o = zeros(nGrid,nGrid);
            for k = 1:param.n_action
                a = param.a_space(k);
                for i = 1 : nGrid
                    for j = 1: nGrid
                        if j == 1
                            f_o(i,j) = normcdf((o_grid(j) + w_omega / 2 - gam_omega0 - ...
                                gam_omega1 * o_grid(i) - gamma_a * a )/sqrt(sigma_omega));
                        elseif j == nGrid
                            f_o(i,j) = 1 - normcdf((o_grid(j) - w_omega / 2 - gam_omega0...
                                - gam_omega1 * o_grid(i) - gamma_a * a)/sqrt(sigma_omega));
                        else
                            f_o(i,j) = normcdf((o_grid(j) + w_omega / 2 - gam_omega0 - ...
                                gam_omega1 * o_grid(i) - gamma_a * a)/sqrt(sigma_omega)) - ...
                                normcdf((o_grid(j) - w_omega / 2 - gam_omega0 - ...
                                gam_omega1 * o_grid(i) - gamma_a * a)/sqrt(sigma_omega));
                        end
                    end
                end
                P{k} = kron(F_z,f_o);
            end
            [tmp0, tmp1,tmp2,tmp3,tmp4] = ndgrid(o_grid,z_grid,z_grid,z_grid,z_grid);
            state = [tmp4(:),tmp3(:),tmp2(:),tmp1(:),tmp0(:)];
        end %end of state transition

        function [at,yt,zt] = simdata(theta_vec, S,nT,nM)
            %---------------------------------------------------------------------    
            % SYNTAX: [at,yt,zt,zt_1] = simdata(theta_vec, S,nT,nM)
            % INPUT: param
            %      P:   transition of exogenous data
            %      MC:  number of montecarlo
            %---------------------------------------------------------------------
            ev = zeros(S.n_state,S.n_action);
            pi = DDCMixture.dpidth(S) * theta_vec;
            [p1,ev] = DDCMixture.solveNFXP(ev,pi,S);
            n_state = S.n_state; %number of exogenous states
            p1 = reshape(p1,n_state,2);
            tmp11= S.P{1} .* repmat(p1(:,1),1,n_state); %y = 1, a = 1
            tmp10= S.P{2} .* repmat((1 - p1(:,1)),1,n_state); %y = 1, a = 0
            tmp01= S.P{1} .* repmat(p1(:,2),1,n_state);%y = 0, a = 1
            tmp00= S.P{2} .* repmat((1 - p1(:,2)),1,n_state);%y = 0, a = 0
            F_x  = [tmp11 , tmp10 ; tmp01, tmp00]; %joint transition of (y,z)
            % ---------------------------------------------------------------------
            %** Data structure
            %        at:     record  a_t for t=1:nT
            %        yt:     record  y_{t-1} for t=1:nT
            %        zt:     record  z_t for t=1:nT
            %        zt_1:   record  z_{t-1} for t=1:nT
            %---------------------------------------------------------------------
            at = zeros(nM, nT);
            yt = zeros(nM, nT);
            zt = zeros(nM, nT);
            zt_1 = zeros(nM, nT); %z t-1
            % ---------------------------------------------------------------------
            %** Create Markov Process of Exogeneous Variable
            % x_t = (y_t,z_t) 
            % need to simulate nT + 2 times
            % need a_t : t = 0 to nT
            %            y_t = a_t-1
            % need z_t : t = 0 to nT
            %---------------------------------------------------------------------
            x_mc  = dtmc(F_x); 
            
            for m = 1:nM
                   x_sim = simulate(x_mc,nT + 1);    
                   a     = floor((x_sim - 1) / n_state) + 1;
                   z     = rem(x_sim - 1, n_state)+1;
                   at(m,:) = a(2:nT+1);
                   yt(m,:) = a(1:nT);
                   zt(m,:) = z(1:nT);
                   
            end
        end
        
        function [at,yt,zt,zt_1] = simdata_mix(S1,S2,theta_vec,param,mix_param)
            %---------------------------------------------------------------------    
            % SYNTAX: [at,yt,zt,zt_1] = simdata(P,state,p1, param)
            % INPUT: param
            %      P:   transition of exogenous data
            %      MC:  number of montecarlo
            %---------------------------------------------------------------------
            %Solve for type I and type II CCP
            %Simulate data
            pd = makedist('Multinomial','probabilities',mix_param.mixProb);
            ubs_type = random(pd,param.nM,1);
            mix_prob = tabulate(ubs_type);
            mix_prob = mix_prob(1,3)/100;
            
            %---------------------------------------------------------------------
            %** Set parameters
            %---------------------------------------------------------------------
            
            [at1,yt1,zt1] = DDCMixture.simdata(theta_vec,S1,param.nT,round(param.nM * mix_prob));
            [at2,yt2,zt2] = DDCMixture.simdata(theta_vec,S2,param.nT,round(param.nM * (1 - mix_prob)));
            at = [at1;at2];
            yt = [yt1;yt2];
            zt = [zt1;zt2];
            
        end %simdata_mix
        
        %---------------------------------------------------------------------    
      
        function [ll,ev,score] = ll(datasim,theta_vec,S)
            %[ll,ev,score] = DDCMixture.ll(datasim,theta_vec,S)
            
            % Log likelihood in NFXP
            % First solve for the equilibrium EV
            a = reshape(datasim.at,S.N*S.T,1);
            z = reshape(datasim.zt,S.N*S.T,1);
            y = reshape(datasim.yt,S.N*S.T,1);
            [ev,p1,F] = NFXP.solve(theta_vec,S);
            n_state = S.n_state;
            ind = (y - 1) * n_state + z ;
            data_p1 = p1(ind);

            ll  = (a==1).*log(data_p1) + (a==2).*log(1 - data_p1);
%             ll  = reshape(ll,S.N,S.T);
            %MULTIPLY THE WEIGHT HERE
%             ll  = reshape(ll,S.N*S.T,1);
            
            [devdth,dpidth] = DDCMixture.devdth(S,F,p1);
            % Step 4: Compute the derivative of log-likelihood w.r.t theta
            score    = bsxfun(@times,( (a == 1) - p1(ind)) , dpidth(ind,:) +...
                      S.beta *  devdth(z,:));  
        end
        
        %---------------------------------------------------------------------    
        
        function [devdth,dpidth] = devdth(S,F,p1)
            % Compute score function
                n_state = S.n_state;
                ccp = reshape(p1,n_state,2);
                % Step 1: Compute the derivative of pi w.r.t to mp
                dvpdth   = exp(S.state(:,5)).* [ones(n_state,1),S.state(:,1),S.state(:,2)]; %Derivative of variable profit
                dfcdth   = [ones(n_state,1),S.state(:,3)];% Derivative of fixed cost
                decdth   = [ones(n_state,1),S.state(:,4)];
                % Step 2: Compute the derivative of contraction mapping w.r.t mp
                dpidth_1 =  [dvpdth,-dfcdth,zeros(size(decdth))];
                dpidth_0 =  [dvpdth,-dfcdth,-decdth];
                dtdth_1  =  ccp(:,1) .* dpidth_1 ; %dpi/dtheta when y = 1
                dtdth_0  =  ccp(:,2) .* dpidth_0 ; %dpi/dtheta when y = 0
                dtdth    = [dtdth_1;dtdth_0];
                dpidth   = [dpidth_1;dpidth_0];

                % This is different than what I've seen, But I don't care 
                % Step 3: Compute the derivative of EV w.r.t theta
                dvdth   = F \ dtdth; 
                dvdth_1 = dvdth(1:n_state,:);
                dvdth_0 = dvdth((1+n_state):(2*n_state),:);
                devdth  = S.P * ( dvdth_1 - dvdth_0 );
        end
        
        function dpidth = dpidth(S)
            n_state = S.n_state;
                % Step 1: Compute the derivative of pi w.r.t to mp
                dvpdth   = exp(S.state(:,5)).* [ones(n_state,1),S.state(:,1),S.state(:,2)]; %Derivative of variable profit
                dfcdth   = [ones(n_state,1),S.state(:,3)];% Derivative of fixed cost
                decdth   = [ones(n_state,1),S.state(:,4)];
                % Step 2: Compute the derivative of contraction mapping w.r.t mp
                dpidth_1 =  [dvpdth,-dfcdth,zeros(size(decdth))];
                dpidth_0 =  [dvpdth,-dfcdth,-decdth];
                dpidth   = [dpidth_1;dpidth_0];
                
        end
    
        
        %------------------------------------------------------------------
        %*************      ESTIMATE DATA PROB         ********************
        %------------------------------------------------------------------
        
        %------------------------------------------------------------------
        %*************         GMM Estimator           ********************
        %------------------------------------------------------------------
        function p1 = gmm(datasim,S)
            eul = 0.5772156649;
            %*************      ESTIMATE DATA PROB         ********************
            n_state = S.n_state;
            a = datasim.at(:);
            z = datasim.zt(:);
            y = datasim.yt(:);
            ind = (y - 1) * n_state + z ;
            freq_tab = crosstab(a,ind);
            p1 = min(max(freq_tab(1,:) ./ (sum(freq_tab,1)),1e-8),1- 1e-8);
            p1 = p1';
%             p1 = reshape(p1,n_state,2);
            
            F_size = size(S.P{2});
            F1 = [S.P{1},zeros(F_size);S.P{1},zeros(F_size)];
            F0 = [zeros(F_size),S.P{2};zeros(F_size),S.P{2}];
            F_opt = (eye(F_size(1)*2) - S.beta * F1) * inv(eye(F_size(1)*2) - S.beta * F0);
            RHS = F_opt * (eul - log(1 - p1(:))) + log(p1(:)) - eul ;
            LHS = DDCMixture.dpidth(S);
            X   = LHS(ind,:);
            Y   = RHS(ind);
            theta_hat = (X'*X)\(X'*Y);
        end
        
        %------------------------------------------------------------------
        %*************      LIKELIHOOD FUNCTIONS       ********************
        %------------------------------------------------------------------
        function [logl, p1,p2,ll_1,ll_2 ] =  ll_NFXP(theta_vec,datasim,w,p1,p2,S1,S2)
            % Likelihood using nested fixed point
            a = reshape(datasim.at,S1.N*S1.T,1);
            z = reshape(datasim.zt,S1.N*S1.T,1);
            y = reshape(datasim.yt,S1.N*S1.T,1);
            [n_state,d_state] = size(S1.state);
            ind = (y - 1) * n_state + z ;
            
            pi_1       = DDCMixture.dpidth(S1) * theta_vec ;
            pi_2       = DDCMixture.dpidth(S2) * theta_vec ;
            
            ev_1 = zeros(S1.n_state,S1.n_action);
            ev_2 = zeros(S2.n_state,S2.n_action);
            [p1_1,ev_1] = DDCMixture.solveNFXP(ev_1,pi_1,S1);
            p1_1 = p1_1(:);
            data_p_1 = p1_1(ind);
            
            [p1_2,ev_2] = DDCMixture.solveNFXP(ev_2,pi_2,S2);
            p1_2 = p1_2(:);
            data_p_2 = p1_2(ind);
               
            ll_1  = sum(reshape([a==1].*log(data_p_1) + [a==2].*log(1 - data_p_1),S1.N,S1.T),2);
            ll_2  = sum(reshape([a==1].*log(data_p_2) + [a==2].*log(1 - data_p_2),S1.N,S1.T),2);
            
            lik = w(:,1) .* ll_1 + w(:,2) .* ll_2 ;
            logl = sum(lik);

            if nargout > 1
                p1 = reshape(p1_1,n_state * 2,1);
                p2 = reshape(p1_2,n_state * 2,1);
            end
        end
        
        function [logl, p1_1,p2_1,ll_1,ll_2 ] =  ll_EE(theta_vec,datasim,w,p1_1,p2_1,S1,S2,invF)
            eul = 0.5772156649;
            a = reshape(datasim.at,S1.N*S1.T,1);
            z = reshape(datasim.zt,S1.N*S1.T,1);
            y = reshape(datasim.yt,S1.N*S1.T,1);
            [n_state,d_state] = size(S1.state);
            ind = (y - 1) * n_state + z ;
            
            pi_1       = DDCMixture.dpidth(S1) * theta_vec ;
            pi_2       = DDCMixture.dpidth(S2) * theta_vec ;
            v_til_1 = pi_1 + S1.beta * [S1.P{1},-S1.P{2};S1.P{1},-S1.P{2}] * invF * (eul - log(1 - p1_1(:) ));
            v_til_2 = pi_2 + S2.beta * [S2.P{1},-S2.P{2};S2.P{1},-S2.P{2}] * invF * (eul - log(1 - p2_1(:) ));
            p1_1 = exp(v_til_1) ./ (1 + exp(v_til_1));
            p2_1 = exp(v_til_2) ./ (1 + exp(v_til_2));
            data_p_1 = p1_1(ind);
            data_p_2 = p2_1(ind);
            ll_1  = sum(reshape([a==1].*log(data_p_1) + [a==2].*log(1 - data_p_1),S1.N,S1.T),2);
            ll_2  = sum(reshape([a==1].*log(data_p_2) + [a==2].*log(1 - data_p_2),S1.N,S1.T),2);
            
            lik = w(:,1) .* ll_1 + w(:,2) .* ll_2 ;
            logl = sum(lik);
        end
         
        function [p1,ev] = mapping_seq(ev,pi,S)
            % Calculate ccp and ev based on sequential mapping
            beta = S.beta;
            P = S.P;
            ev = reshape(ev,size(P{1},2),2);
            pi = reshape(pi,size(P{1},2),2);
            eul = 0.5772156649;
            v_til(:,1) =  pi(:,1) + beta *(P{1}*ev(:,1) - P{2}*ev(:,2)) ;
            v_til(:,2) =  pi(:,2) + beta *(P{1}*ev(:,1) - P{2}*ev(:,2));
            v_til = v_til(:);
            p1 = (exp(v_til) ./ (1 + exp(v_til)));
%             tmp(:,1) = max( pi(:,1) + beta * P *  ev(:,1),  beta * P *  ev(:,2)) ;
%             tmp(:,2) = max( pi(:,2) + beta * P *  ev(:,1),  beta * P *  ev(:,2)) ;            
            tmp(:,1) = eul + log( exp( pi(:,1) + beta * P{1} * ( ev(:,1))) + exp( beta * P{2} *  ev(:,2))) ;
            tmp(:,2) = eul + log( exp( pi(:,2) + beta * P{1} * ( ev(:,1))) + exp( beta * P{2} *  ev(:,2))) ;
            ev = tmp;
        end
        
        function [p1,ev] = mapping_seq_ee(ev,pi,S)
            % Calculate ccp and ev based on sequential mapping
            beta = S.beta;
            P = S.P;
            ev = reshape(ev,size(P{1},2),2);
            pi = reshape(pi,size(P{1},2),2);
            eul = 0.5772156649;
            v_til(:,1) =  pi(:,1) + beta *(P{1}*ev(:,1) - P{2}*ev(:,2));
            v_til(:,2) =  pi(:,2) + beta *(P{1}*ev(:,1) - P{2}*ev(:,2));
            p1 = (exp(v_til) ./ (1 + exp(v_til)));       
            tmp(:,1) = eul -  log( 1 - p1(:,1))  +  beta * P{2} *  ev(:,2) ;
            tmp(:,2) = eul -  log( 1 - p1(:,2))  +  beta * P{2} *  ev(:,2) ;
            p1 = p1(:);
            ev = tmp;
        end
        
        function [p1,ev] = solveNFXP(ev,pi,S)
            tol = 1e-8;
            diff = inf;
            max_iter = 100;
            iter = 0;
            while (diff > tol) 
                [p1_1, ev_1] = DDCMixture.mapping_seq(ev,pi,S);
                diff = max(abs(ev_1 - ev));
                ev = ev_1;
                p1 = p1_1;
                iter = iter + 1;
                if (iter > max_iter)
                    break;
                end
            end
        end
        
        function [logl,ev_1,ev_2,ll_1,ll_2] = ll_SEQ(theta_vec,datasim,w,ev_1,ev_2,S1,S2)
            % This function gives sum of log-likelihood when fixing the p1(h-1)
            a = reshape(datasim.at,S1.N*S1.T,1);
            z = reshape(datasim.zt,S1.N*S1.T,1);
            y = reshape(datasim.yt,S1.N*S1.T,1);
            [n_state,d_state] = size(S1.state);
            ind = (y - 1) * n_state + z ;
            
            pi_1       = DDCMixture.dpidth(S1) * theta_vec ;
            pi_2       = DDCMixture.dpidth(S2) * theta_vec ;
            
            [p1_1,ev_1] =  DDCMixture.mapping_seq(ev_1,pi_1,S1);
            p1_1 = p1_1(:);
            data_p_1 = p1_1(ind);
            
            [p1_2,ev_2] =  DDCMixture.mapping_seq(ev_2,pi_2,S2);
            p1_2 = p1_2(:);
            data_p_2 = p1_2(ind);
               
            ll_1  = sum(reshape([a==1].*log(data_p_1) + [a==2].*log(1 - data_p_1),S1.N,S1.T),2);
            ll_2  = sum(reshape([a==1].*log(data_p_2) + [a==2].*log(1 - data_p_2),S1.N,S1.T),2);
            

            lik = w(:,1) .* ll_1 + w(:,2) .* ll_2 ;
            logl = sum(lik);

            if nargout > 1
                p1 = reshape(p1_1,n_state * 2,1);
                p2 = reshape(p1_2,n_state * 2,1);
            end
        end

        function [logl,ev_1,ev_2,ll_1,ll_2] = ll_SEQ_EE(theta_vec,datasim,w,ev_1,ev_2,S1,S2)
            % This function gives sum of log-likelihood when fixing the p1(h-1)
            a = reshape(datasim.at,S1.N*S1.T,1);
            z = reshape(datasim.zt,S1.N*S1.T,1);
            y = reshape(datasim.yt,S1.N*S1.T,1);
            [n_state,d_state] = size(S1.state);
            ind = (y - 1) * n_state + z ;
            
            pi_1       = DDCMixture.dpidth(S1) * theta_vec ;
            pi_2       = DDCMixture.dpidth(S2) * theta_vec ;
            
            [p1_1,ev_1] =  DDCMixture.mapping_seq_ee(ev_1,pi_1,S1);
            p1_1 = p1_1(:);
            data_p_1 = p1_1(ind);
            
            [p1_2,ev_2] =  DDCMixture.mapping_seq_ee(ev_2,pi_2,S2);
            p1_2 = p1_2(:);
            data_p_2 = p1_2(ind);
               
            ll_1  = sum(reshape([a==1].*log(data_p_1) + [a==2].*log(1 - data_p_1),S1.N,S1.T),2);
            ll_2  = sum(reshape([a==1].*log(data_p_2) + [a==2].*log(1 - data_p_2),S1.N,S1.T),2);
            
            lik = w(:,1) .* ll_1 + w(:,2) .* ll_2 ;
            logl = sum(lik);

            if nargout > 1
                p1 = reshape(p1_1,n_state * 2,1);
                p2 = reshape(p1_2,n_state * 2,1);
            end
        end
        
        function [logl, p1,p2,ll_1,ll_2 ] =  ll_FD(theta_vec,datasim,w,p1,p2,S1,S2)
            %update using finite dependence
            a = reshape(datasim.at,S1.N*S1.T,1);
            z = reshape(datasim.zt,S1.N*S1.T,1);
            y = reshape(datasim.yt,S1.N*S1.T,1);
            [n_state,d_state] = size(S1.state);
            ind = (y - 1) * n_state + z ;
            ccp_1 = reshape(p1,n_state,2);
            ccp_2 = reshape(p2,n_state,2);
            
            pi_1       = DDCMixture.dpidth(S1) * theta_vec ;
            pi_2       = DDCMixture.dpidth(S2) * theta_vec ;
            
            tmp_1 = pi_1 + S1.beta  * kron([1;1],(S1.P{1} * log(1- ccp_1(:,2)) - S1.P{2} * log(1 - ccp_1(:,1)) ));
            p1_1 = exp(tmp_1) ./ ( 1 + exp(tmp_1));
            data_p_1 = p1_1(ind);
            
            tmp_2 = pi_2 + S2.beta * kron([1;1],  (S2.P{1} * log(1- ccp_2(:,2)) - S2.P{1} *log(1 - ccp_2(:,1)) ));
            p1_2 = exp(tmp_2) ./ ( 1 + exp(tmp_2));
            data_p_2 = p1_2(ind);
               
            ll_1  = sum(reshape([a==1].*log(data_p_1) + [a==2].*log(1 - data_p_1),S1.N,S1.T),2);
            ll_2  = sum(reshape([a==1].*log(data_p_2) + [a==2].*log(1 - data_p_2),S1.N,S1.T),2);
            

            lik = w(:,1) .* ll_1 + w(:,2) .* ll_2 ;
            logl = sum(lik);

            if nargout > 1
                p1 = reshape(p1_1,n_state * 2,1);
                p2 = reshape(p1_2,n_state * 2,1);
            end
        end
        
        
        
        %------------------------------------------------------------------
        % **************      SEQUENTIAL ESTIMATION       *****************
        %------------------------------------------------------------------

        function [theta_hat,w,iter] = SequentialEstimation(datasim,mix_param,S1,S2,theta_vec0,opt)
            o1 = optimset('LargeScale','off','Display','off');
            o2 = optimset('LargeScale','off','Display','off','GradObj','on');

            %STEP 1: guess the weight for the two type mixture
            w0 = [0.1 * ones(S1.N,1), 0.9 * ones(S1.N,1) ]; %The weight
            q = mean(w0); %The probability of each type
            
            opt.tol  = 1e-3;
            opt.max_iter = 100;
            opt.output = 0;
            diff      = 1;
            iter      = 0;
            if string(opt.method) == 'FD'
%                 fprintf('Solving using Finite Dependence\n');
                p1_1 = ones(S1.n_state,S1.n_action)/S1.n_action;
                p2_1 = ones(S2.n_state,S2.n_action)/S2.n_action;
                ll_function   = @DDCMixture.ll_FD;
                while (diff > opt.tol)
                    ts = tic;
                    %STEP 2: Then estimate the theta using the mixture weight
                    f = @(theta_vec)(-ll_function(theta_vec,datasim,w0,p1_1,p2_1,S1,S2)); %Make sure the solver works
                    theta_vec = fmincon(f,theta_vec0,[],[],[],[],[],[],[],o1);
                    if opt.output > 0
                        fprintf('Iteration: %d, Difference: %.4f, Time elapsed: %f Seconds \n',iter,diff,Result.REALtime);
                    end                %STEP3: Update the weight using the theta
                    %STEP 3: Update weight and q
                    [ll,p1_1,p2_1,ll_1,ll_2] = ll_function(theta_vec,datasim,w0,p1_1,p2_1,S1,S2);

                    minl = min(min(ll_1,ll_2));

                    w1(:,1) = (q(1) * exp(ll_1 - minl) );
                    w1(:,2) = (q(2) * exp(ll_2 - minl) );
                    w1 = w1 ./sum(w1,2);
                    diff = max(max(abs(w0 - w1)));
                    w0 = w1;
                    q = mean(w1); %The probability of each type
                    % STEP 4: Update conditional choice probability
                    theta_vec0 = theta_vec;
                    iter = iter + 1;
                    if iter > opt.max_iter
                        break;
                    end
                end
            elseif string(opt.method) == 'EE'
                F_size = size(S1.P{2});
                F0 = [zeros(F_size),S1.P{2};zeros(F_size),S1.P{2}];
                size_n_state = S1.n_state * S1.n_action;
                invF =  inv(eye(size_n_state) - S1.beta * F0);
                p1_1 = ones(S1.n_state,S1.n_action)/S1.n_action;
                p2_1 = ones(S2.n_state,S2.n_action)/S2.n_action;
                ll_function   = @DDCMixture.ll_EE;
%                 fprintf('Solving using Euler Equation\n');

                while (diff > opt.tol)
                    ts = tic;
                    f = @(theta_vec)(-ll_function(theta_vec,datasim,w0,p1_1,p2_1,S1,S2,invF)); 
                    theta_vec = fmincon(f,theta_vec0,[],[],[],[],[],[],[],o1);
                    if opt.output > 0
                        fprintf('Iteration: %d, Difference: %.4f, Time elapsed: %f Seconds \n',iter,diff,Result.REALtime);
                    end
                    [ll,p1_1,p2_1,ll_1,ll_2] = ll_function(theta_vec,datasim,w0,p1_1,p2_1,S1,S2,invF);

                    minl = min(min(ll_1,ll_2));

                    w1(:,1) = (q(1) * exp(ll_1 - minl) );
                    w1(:,2) = (q(2) * exp(ll_2 - minl) );
                    w1 = w1 ./sum(w1,2);
                    diff = max(max(abs(w0 - w1)));
                    w0 = w1;
                    q = mean(w1); %The probability of each type
                    % STEP 4: Update conditional choice probability
                    theta_vec0 = theta_vec;
                    iter = iter + 1;
                    if iter > opt.max_iter
                        break;
                    end
                end

                
            elseif string(opt.method) == 'NFXP'
                %Solve initial CCP
                ev_1 = zeros(S1.n_state,S1.n_action);
                ev_2 = zeros(S2.n_state,S2.n_action);
                p1_1 = ones(S1.n_state,S1.n_action)/S1.n_action;
                p2_1 = ones(S2.n_state,S2.n_action)/S2.n_action;
%                 fprintf('Solving using Nested fixed point method\n');
                ll_function   = @DDCMixture.ll_NFXP;
                
                while (diff > opt.tol)
                    ts = tic;
                    f = @(theta_vec)(-ll_function(theta_vec,datasim,w0,p1_1,p2_1,S1,S2)); 
                    theta_vec = fmincon(f,theta_vec0,[],[],[],[],[],[],[],o1);
                    if opt.output > 0
                        fprintf('Iteration: %d, Difference: %.4f, Time elapsed: %f Seconds \n',iter,diff,Result.REALtime);
                    end
                    [ll,p1_1,p2_1,ll_1,ll_2] = ll_function(theta_vec,datasim,w0,p1_1,p2_1,S1,S2);

                    minl = min(min(ll_1,ll_2));

                    w1(:,1) = (q(1) * exp(ll_1 - minl) );
                    w1(:,2) = (q(2) * exp(ll_2 - minl) );
                    w1 = w1 ./sum(w1,2);
                    diff = max(max(abs(w0 - w1)));
                    w0 = w1;
                    q = mean(w1); %The probability of each type
                    % STEP 4: Update conditional choice probability
                    theta_vec0 = theta_vec;
                    iter = iter + 1;
                    if iter > opt.max_iter
                        break;
                    end
                end
            elseif string(opt.method) == 'SEQ'                 
                ev_1 = zeros(S1.n_state,S1.n_action);
                ev_2 = zeros(S2.n_state,S2.n_action);
                ll_function = @DDCMixture.ll_SEQ;
%                 fprintf('Solving using Sequential Estimation\n');
                
                while (diff > opt.tol)
                    ts = tic;
                    f = @(theta_vec)(-ll_function(theta_vec,datasim,w0,ev_1,ev_2,S1,S2)); %Make sure the solver works
                    theta_vec = fmincon(f,theta_vec0,[],[],[],[],[],[],[],o1);
                    if opt.output > 0
                        fprintf('Iteration: %d, Difference: %.4f, Time elapsed: %f Seconds \n',iter,diff,Result.REALtime);
                    end
                    [ll,ev_1,ev_2,ll_1,ll_2] = ll_function(theta_vec,datasim,w0,ev_1,ev_2,S1,S2);
                    
                    minl = min(min(ll_1,ll_2));

                    w1(:,1) = (q(1) * exp(ll_1 - minl) );
                    w1(:,2) = (q(2) * exp(ll_2 - minl) );
                    w1 = w1 ./sum(w1,2);
                    diff = max(max(abs(w0 - w1)));
                    w0 = w1;
                    q = mean(w1); %The probability of each type
                    % STEP 4: Update conditional choice probability
                    theta_vec0 = theta_vec;
                    
                    iter = iter + 1;
                    if iter > opt.max_iter
                        break;
                    end
                end
                elseif string(opt.method) == 'SEQ_EE'                 
                ev_1 = zeros(S1.n_state,S1.n_action);
                ev_2 = zeros(S2.n_state,S2.n_action);
                ll_function = @DDCMixture.ll_SEQ_EE;
%                 fprintf('Solving using Sequential Estimation\n');
                
                while (diff > opt.tol)
                    ts = tic;
                    f = @(theta_vec)(-ll_function(theta_vec,datasim,w0,ev_1,ev_2,S1,S2)); %Make sure the solver works
                    theta_vec = fmincon(f,theta_vec0,[],[],[],[],[],[],[],o1);
                    if opt.output > 0
                        fprintf('Iteration: %d, Difference: %.4f, Time elapsed: %f Seconds \n',iter,diff,Result.REALtime);
                    end
                    [ll,ev_1,ev_2,ll_1,ll_2] = ll_function(theta_vec,datasim,w0,ev_1,ev_2,S1,S2);
                    
                    minl = min(min(ll_1,ll_2));

                    w1(:,1) = (q(1) * exp(ll_1 - minl) );
                    w1(:,2) = (q(2) * exp(ll_2 - minl) );
                    w1 = w1 ./sum(w1,2);
                    diff = max(max(abs(w0 - w1)));
                    w0 = w1;
                    q = mean(w1); %The probability of each type
                    % STEP 4: Update conditional choice probability
                    theta_vec0 = theta_vec;
                    
                    iter = iter + 1;
                    if iter > opt.max_iter
                        break;
                    end
                end
            else
                fprintf('No correct algorithm chosen\n');
            end
            w = w1;
            theta_hat = theta_vec;
        end %End of Mixture.sequential_EE

        
        %------------------------------------------------------------------
        % **************                                  *****************
        % **************      Single Type Estimation      *****************
        % **************                                  *****************
        %------------------------------------------------------------------
        function [logl, p1,ll_1]  =  ll_EE_s(theta_vec,datasim,p1,S1,invF)
            eul = 0.5772156649;
            a = reshape(datasim.at,S1.N*S1.T,1);
            z = reshape(datasim.zt,S1.N*S1.T,1);
            y = reshape(datasim.yt,S1.N*S1.T,1);
            [n_state,d_state] = size(S1.state);
            ind = (y - 1) * n_state + z ;
            F_0 = kron([0,1;0,1],S1.P{2});
            F_1 = kron([1,0;1,0],S1.P{1});
            F_til = F_1 - F_0;
            
            p1 = vec(p1);
            pi_1       = DDCMixture.dpidth(S1) * theta_vec ;
            
            % v_til_1 = pi_1 + S1.beta * [S1.P{1},-S1.P{2};S1.P{1},-S1.P{2}] * invF * (eul - log(1 - p1_1(:) ));
            tmp_1 = pi_1 + S1.beta  * F_til *  invF * (eul - log(1 - p1)); 
            
            p1_1 = exp(tmp_1) ./ (1 + exp(tmp_1));
            data_p_1 = p1_1(ind);
            ll_1  = sum(reshape([a==1].*log(data_p_1) + [a==2].*log(1 - data_p_1),S1.N,S1.T),2);
           
            logl = sum(ll_1);
            if nargout > 1
                p1 = reshape(p1_1,n_state * 2,1);
            end
        end
        
        function [logl, p1,ll_1]  =  ll_HM_s(theta_vec,datasim,p1,S1)
            eul = 0.5772156649;
            a = reshape(datasim.at,S1.N*S1.T,1);
            z = reshape(datasim.zt,S1.N*S1.T,1);
            y = reshape(datasim.yt,S1.N*S1.T,1);
            [n_state,d_state] = size(S1.state);
            ind = (y - 1) * n_state + z ;
            F_0 = kron([0,1;0,1],S1.P{2});
            F_1 = kron([1,0;1,0],S1.P{1});
            F_til = F_1 - F_0;
            
            p1 = vec(p1);
            pi_1       = DDCMixture.dpidth(S1) * theta_vec ;
            pi_P       = eul + ( pi_1 -  log(p1) ).*  p1 - log(1 - p1 ) .*(1- p1);
            F_P        = diag(p1) * F_til + F_0;
            invF       = inv(eye(size(F_P,2)) - S1.beta * F_P);
            % v_til_1 = pi_1 + S1.beta * [S1.P{1},-S1.P{2};S1.P{1},-S1.P{2}] * invF * (eul - log(1 - p1_1(:) ));
            tmp_1 = pi_1 + S1.beta  * F_til *  invF *  pi_P ; 
            
            p1_1 = exp(tmp_1) ./ (1 + exp(tmp_1));
            data_p_1 = p1_1(ind);
            ll_1  = sum(reshape([a==1].*log(data_p_1) + [a==2].*log(1 - data_p_1),S1.N,S1.T),2);
           
            logl = sum(ll_1);
            if nargout > 1
                p1 = reshape(p1_1,n_state * 2,1);
            end
        end
        
        function [logl, p1,ll_1] =  ll_FD_s(theta_vec,datasim,p1,S1,p_star,V_star)
            %update using finite dependence
            if ~exist("V_star"); V_star = zeros(S1.n_state * S1.n_action,1); end;
            a = reshape(datasim.at,S1.N*S1.T,1);
            z = reshape(datasim.zt,S1.N*S1.T,1);
            y = reshape(datasim.yt,S1.N*S1.T,1);
            [n_state,d_state] = size(S1.state);
            ind = (y - 1) * n_state + z ;
            
            F_0 = kron([0,1;0,1],S1.P{2});
            F_1 = kron([1,0;1,0],S1.P{1});
            F_til = F_1 - F_0;
%             F   = [F_0;F_1];
            
            p1 = vec(p1);
            pi_1       = DDCMixture.dpidth(S1) * theta_vec ;
            F_star = diag(p_star) * F_til + F_0;

            % tmp_1 = pi_1 + S1.beta  * kron([1;1],(S1.P{1} * log(1- ccp_1(:,2)) - S1.P{2} * log(1 - ccp_1(:,1)) ));
            tmp_1 = pi_1 + S1.beta  * F_til * ( p_star .* pi_1 +   0.577215665 - (1 - p_star) .* log(1 - p1) - p_star .* log(p1) + S1.beta * F_star * V_star);
            
            p1_1 = exp(tmp_1) ./ ( 1 + exp(tmp_1));
            
            data_p_1 = p1_1(ind);
            
               
            ll_1  = sum(reshape([a==1].*log(data_p_1) + [a==2].*log(1 - data_p_1),S1.N,S1.T),2);           

            lik = ll_1 ;
            logl = sum(lik);

            if nargout > 1
                p1 = reshape(p1_1,n_state * 2,1);
            end
        end
        
        % **************      Function *****************
        function p1 = estimate_ccp(datasim,S1)
            %p1 = [p(1|z,1),p(1|z,0)]
            a = datasim.at(:);
            z = datasim.zt(:);
            y = datasim.yt(:);
            ind = (y - 1) *S1.n_state + z;
            
            for each_ind = 1:max(ind)
                subset_ind = (ind == each_ind);
                if sum(subset_ind) > 0
                    p1(each_ind) = mean(a(subset_ind)==1);
                    if p1(each_ind)==0
                        p1(each_ind) = 1/length(a);
                    elseif p1(each_ind)==1
                        p1(each_ind) = 1-1/length(a);
                    end
                else 
                    p1(each_ind) = 1/length(a);
                end
            end
            p1 = vec(p1);
        end
        
        
        function V = bellman_ee(V,pi_1,p1_1,F1,F0,beta)
         	V = p1_1 .* (pi_1 + log(p1_1) + beta  * F1 * V) + (1-p1_1).* (log(1-p1_1) + beta  * F0 * V);
        end
        
        
        function [theta_hat,iter] = SingleEstimation(datasim,S1,theta_vec0,p_star,opt)
            o1 = optimset('LargeScale','off','Display','off');
            o2 = optimset('LargeScale','off','Display','off','GradObj','on');

            
            opt.tol  = 1e-8;
%             opt.max_iter = 100;
            opt.output = 0;
            diff      = 1;
            iter      = 0;
            if ~isfield("true_ccp",opt)
                opt.true_ccp=0;
            end
            if opt.true_ccp
                p1_1 = S1.p_1;
            else
                p1_1 = DDCMixture.estimate_ccp(datasim,S1); %The initial CCP
            end
            
            if string(opt.method) == 'FD'
                
                ll_function   = @DDCMixture.ll_FD_s;

                while (diff > opt.tol)&(iter<opt.max_iter)
                    
                    %STEP 2: Then estimate the theta using the mixture weight
                    f = @(theta_vec)(-ll_function(theta_vec,datasim,p1_1,S1,p_star)); %Make sure the solver works

                    theta_vec = fmincon(f,theta_vec0,[],[],[],[],[],[],[],o1);
                    if opt.output > 0
                        fprintf('Iteration: %d, Difference: %.4f, Time elapsed: %f Seconds \n',iter,diff,Result.REALtime);
                    end                %STEP3: Update the weight using the theta
                    %STEP 3: Update weight and q
                    [ll,p1,ll_1] = ll_function(theta_vec,datasim,p1_1,S1,p_star);

                    diff = max(abs(vec(p1_1) - vec(p1) ));

                    % STEP 4: Update conditional choice probability
                    theta_vec0 = theta_vec;
                    p1_1 = p1;
                    iter = iter + 1;
                    if iter > opt.max_iter
                        break;
                    end
                end
                
            elseif string(opt.method) == 'AFD2'
                
                ll_function   = @DDCMixture.ll_FD_s;
                F_size = size(S1.P{2});
                    
                %STEP 2: Then estimate the theta using the mixture weight
                f = @(theta_vec)(-ll_function(theta_vec,datasim,p1_1,S1,p_star)); %Make sure the solver works

                theta_vec = fmincon(f,theta_vec0,[],[],[],[],[],[],[],o1);
                if opt.output > 0
                    fprintf('Iteration: %d, Difference: %.4f, Time elapsed: %f Seconds \n',iter,diff,Result.REALtime);
                end                %STEP3: Update the weight using the theta
                
                % STEP 4: Update conditional choice probability
                [ll,p1,ll_1] = ll_function(theta_vec,datasim,p1_1,S1,p_star);
                theta_vec0 = theta_vec; p1_1 = p1;
                % % V_star = inv(eye(size_n_state) - S1.beta * F0) * (.5772156649 - log(1 - p1)); %use EE to estimate value function
                F0 = [zeros(F_size),S1.P{2};zeros(F_size),S1.P{2}];
                F1 = [S1.P{2},zeros(F_size);S1.P{2},zeros(F_size)];
                size_n_state = S1.n_state * S1.n_action;
                V_star = zeros(size_n_state,1);
                pi_1  = DDCMixture.dpidth(S1) * theta_vec;
                for q = 1:10
                    V_star = DDCMixture.bellman_ee(V_star,pi_1,p1_1,F1,F0,S1.beta);
                end
                
                diff = inf;
                while (diff > opt.tol) &(iter<opt.max_iter)
                    %STEP 2: Then estimate the theta using the mixture weight
                    f = @(theta_vec)(-ll_function(theta_vec,datasim,p1_1,S1,p_star,V_star)); %Make sure the solver works

                    theta_vec = fmincon(f,theta_vec0,[],[],[],[],[],[],[],o1);
                    if opt.output > 0
                        fprintf('Iteration: %d, Difference: %.4f, Time elapsed: %f Seconds \n',iter,diff,Result.REALtime);
                    end
                    [ll,p1,ll_1] = ll_function(theta_vec,datasim,p1_1,S1,p_star);
                    diff = max(abs(vec(p1_1) - vec(p1) ));
                    % STEP 4: Update conditional choice probability
                    
                    pi_1  = DDCMixture.dpidth(S1) * theta_vec;
                    for q = 1:10
                        V_star = DDCMixture.bellman_ee(V_star,pi_1,p1_1,F1,F0,S1.beta);
                    end
                    theta_vec0 = theta_vec; p1_1 = p1;
                    iter = iter + 1;
                end
                
                
            elseif string(opt.method) == 'EE'
                F_size = size(S1.P{2});
                F0 = [zeros(F_size),S1.P{2};zeros(F_size),S1.P{2}];
                size_n_state = S1.n_state * S1.n_action;
                invF =  inv(eye(size_n_state) - S1.beta * F0);
                p1_1 = ones(S1.n_state,S1.n_action)/S1.n_action;
                ll_function   = @DDCMixture.ll_EE_s;
%                 fprintf('Solving using Euler Equation\n');
                
                while (diff > opt.tol) &(iter<opt.max_iter)
                    ts = tic;
                    f = @(theta_vec)(-ll_function(theta_vec,datasim,p1_1,S1,invF)); 
                    theta_vec = fmincon(f,theta_vec0,[],[],[],[],[],[],[],o1);
                    if opt.output > 0
                        fprintf('Iteration: %d, Difference: %.4f, Time elapsed: %f Seconds \n',iter,diff,Result.REALtime);
                    end
                    [ll,p1,ll_1] = ll_function(theta_vec,datasim,p1_1,S1,invF);
                    diff = max(abs(vec(p1_1) - vec(p1) ));
                    % STEP 4: Update conditional choice probability
                    theta_vec0 = theta_vec;
                    p1_1 = p1;
                    iter = iter + 1;
                end
            
            elseif string(opt.method) == 'HM'
                F_size = size(S1.P{2});
                F0 = [zeros(F_size),S1.P{2};zeros(F_size),S1.P{2}];
                size_n_state = S1.n_state * S1.n_action;
                p1_1 = ones(S1.n_state,S1.n_action)/S1.n_action;
                ll_function   = @DDCMixture.ll_HM_s;
%                 fprintf('Solving using Hotz-Miller Inversion\n');
                
                while (diff > opt.tol) &(iter<opt.max_iter)
                    ts = tic;
                    f = @(theta_vec)(-ll_function(theta_vec,datasim,p1_1,S1)); 
                    theta_vec = fmincon(f,theta_vec0,[],[],[],[],[],[],[],o1);
                    if opt.output > 0
                        fprintf('Iteration: %d, Difference: %.4f, Time elapsed: %f Seconds \n',iter,diff,Result.REALtime);
                    end
                    [ll,p1,ll_1] = ll_function(theta_vec,datasim,p1_1,S1);
                    diff = max(abs(vec(p1_1) - vec(p1) ));
                    % STEP 4: Update conditional choice probability
                    theta_vec0 = theta_vec;
                    p1_1 = p1;
                    iter = iter + 1;
                end
                
            end %end of if for type
            theta_hat = theta_vec;

        end %End of Mixture.sequential_EE
    end
end 