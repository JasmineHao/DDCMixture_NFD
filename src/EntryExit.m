% EntryExitGame: procedure to compute contraction operator
%   state : (z_t,y_t) 
%           z = (z_1,z_2,z_3,z_4,omega)
%           
%   action: a_t
%           y_t = a_t-1
% See also: nfxp.bellman, nfxp.dbellman, nfxp.solve, nfxp.maxlik

classdef EntryExit
    methods (Static)
        function Pi = profit(state,y,theta)
            % Input: 
            %       z: vector of exogenous variable
            %          state = [z_1,z_2,z_3,z_4,omega]
            %          Profit Pi =  Variable profit - Fixed cost - Entry cost
            %         
            %       y: endogeneous variable
            %        y = {0,1}
            VP =  (theta.VP0 + theta.VP1 * state(:,1) + theta.VP2 * state(:,2)) .* exp(state(:,5));
            FC = theta.FC0 + theta.FC1 * state(:,3);
            EC = (1 - y) * (theta.EC0 + theta.EC1 * state(:,4));
            Pi = VP - FC - EC;
        end %end of profit
        
        function [P,state] = statetransition(param)
            %--------------------------------------------------------------
            % EntryExit.statetransition
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
                            gam_z1 * z_grid(i))/sqrt(sigma_z));
                    end 
                end
            end
            %--------------------------------------------------------------
            % f_o: transition density for omega_j
            f_o = zeros(nGrid,nGrid);
            for i = 1 : nGrid
                for j = 1: nGrid
                    if j == 1
                        f_o(i,j) = normcdf((o_grid(j) + w_omega / 2 - gam_omega0 - ...
                            gam_omega1 * o_grid(i))/sqrt(sigma_omega));
                    elseif j == nGrid
                        f_o(i,j) = 1 - normcdf((o_grid(j) - w_omega / 2 - gam_omega0...
                            - gam_omega1 * o_grid(i))/sqrt(sigma_omega));
                    else
                        f_o(i,j) = normcdf((o_grid(j) + w_omega / 2 - gam_omega0 - ...
                            gam_omega1 * o_grid(i))/sqrt(sigma_omega)) - ...
                            normcdf((o_grid(j) - w_omega / 2 - gam_omega0 - ...
                            gam_omega1 * o_grid(i))/sqrt(sigma_omega));
                    end 
                end
            end
            F_z = kron(f_z,kron(f_z,kron(f_z,f_z)));
            P = kron(F_z,f_o);
            [tmp0, tmp1,tmp2,tmp3,tmp4] = ndgrid(o_grid,z_grid,z_grid,z_grid,z_grid);
            state = [tmp4(:),tmp3(:),tmp2(:),tmp1(:),tmp0(:)];
        end %end of state transition
        
        function [ev1, p1]=bellman(ev, Pi, P, beta,bparam)
            % NFXP.BELMANN:     Procedure to compute bellman equation
            %
            % Syntax :         ev=EntryExit.bellman(ev, Pi, P, beta,bparam);
            %
            % Inputs :
            %  ev              n_state x 2 matrix of expected values given initial
            %                   guess on value function
            %  Pi              n_state x 2 matrix of one-period payoff
            %  P               n_state x n_state transition matrix of 
            % Output:
            %  ev1              n x 1 matrix of expected values given initial
            %                   guess of ev
            % See also:
            %   nfxp.maxlik, nfxp.solve, nfxp.dbellman
            ev = reshape(ev,bparam.n_state,bparam.n_action);
            ev_11 = Pi(:,1) + beta * P * ev(:,1);
            ev_01 = Pi(:,2) + beta * P * ev(:,1);
            ev_00  = beta * P * ev(:,2);
            ev_1= max(ev_11,ev_00);
            ev_0= max(ev_01,ev_00);
            ev1 = [ev_1;ev_0];
            if nargout>1
                %also compute choice probability from ev (initial input)
                p1 = [1 ./ ( 1 + exp(ev_00 - ev_11)); 1 ./ ( 1 + exp(ev_00 - ev_01))];
            end
        end % end of EntryExit.bellman
        
        
        function dev = dbellman(p1, P, beta)
            % Bellman equation
            % EntryExit.DBELMANN:   Procedure to compute Frechet derivative of Bellman operator
            % Syntax :          ev=nfpx.dbellman(pk,  P, mp);
            n_state = size(P,2);
            p1 = reshape(p1,n_state,2);
            %tmp11 : Prob(a=1|y=1) , tmp10 : Prob(a=0|y=1)
            tmp11= P .* repmat(p1(:,1),1,n_state); %y = 1, a = 1
            tmp10= P .* repmat((1 - p1(:,1)),1,n_state); %y = 1, a = 0
            tmp01= P .* repmat(p1(:,2),1,n_state);%y = 0, a = 1
            tmp00= P .* repmat((1 - p1(:,2)),1,n_state);%y = 0, a = 0
            dev   = beta * [tmp11 , tmp10 ; tmp01, tmp00];
        end 
        
        
        function [ev,p1,F] = solve(P,state,beta,theta)
            %----------------------------------------------------------------------
            %  SYNTAX: [ev,p1] = EntryExit.solve(P,state,beta,theta);
            %----------------------------------------------------------------------
            %  INPUT:
            %     EV0 :      m x 1 matrix or scalar zero. Initial guess of choice specific expected value function, EV.
            %     param:     parameters
            %     state:     State transition matrix
            %     P:         State transition
            %----------------------------------------------------------------------
            %
            %  OUTPUT:
            %     P1:     m x 1 matrix of conditional choice probabilities, P(d=enter|x)
            %     EV:     m x 1 matrix of expected value functions, EV(x)
            %     F:      m x m matrix of derivatives of Identity matrix I minus 
            %             contraction mapping operator, I-T' where T' refers to derivative of the expected  value function
            %----------------------------------------------------------------------
            opt.max_fxpiter= 3;             % Maximum number of times to switch between Newton-Kantorovich iterations and contraction iterations.
            opt.min_cstp   = 4;             % Minimum numbe'r of contraction steps
            opt.max_cstp   = 20;            % Maximum number of contraction steps
            opt.ctol       = 0.02;          % Tolerance before switching to N-K algorithm
            opt.rtol       = 0.02;          % Relative tolerance before switching to N-K algorithm
            opt.nstep      = 20;            % Maximum number of Newton-Kantorovich steps
            opt.ltol0      = 1.0e-10;       % Final exit tolerance in fixed point algorithm, measured in units of numerical precision
            opt.printfxp   = 0;             % print iteration info for fixed point algorithm if > 0
            opt.rtolnk     = .5;             % Tolerance for discarding N-K iteration and move to SA (tolm1/tolm > 1+opt.rtolnk) 
            
            %   n_state: number of states
            %   d_state: dimension of states
            [n_state,d_state] = size(state);
            bparam.n_state = n_state;
            bparam.n_action = 2;
            bparam.d_state = d_state;
            Pi_1 = EntryExit.profit(state,1,theta);
            Pi_0 = EntryExit.profit(state,0,theta);
            Pi   = [Pi_1,Pi_0];
            
            % Initial guess of ev
            ev = zeros(n_state*2,1);
            %Initialize counters
            NKIter=0;
            BellmanIter=0;
            converged=false;%initialize convergence indicator
            tolm=1;
            for k=1:opt.max_fxpiter; %poli-algorithm loop (switching between SA and N-K and back)

                % SECTION A: CONTRACTION ITERATIONS
                if opt.printfxp>0
                    fprintf('\n');
                    fprintf('Begin contraction iterations (for the %d. time)\n',k);
                    fprintf('   j           tol        tol(j)/tol(j-1) \n');
                end;
                %SA contraction steps
                for j=1:opt.max_cstp;
                    ev1=EntryExit.bellman(ev, Pi, P, beta,bparam);

                    BellmanIter=BellmanIter+1;
                    
                    tolm1=max(max(abs(ev-ev1)));
                    rtolm=tolm1/tolm;
                    
                    if opt.printfxp>0
                        fprintf(' %3.0f   %16.8f %16.8f\n',j, tolm1,rtolm);
                    end
                                        %prepare for next iteration
                    ev=ev1;
                    tolm=tolm1;

                    %stopping criteria
                    if (j>=opt.min_cstp) && (tolm1 < opt.ctol)
                        %go to NK iterations due to absolute tolerance
                        break;
                    end;
                    if (j>=opt.min_cstp) && (abs(beta-rtolm) < opt.rtol)
                        %go to NK iterations due to relative tolerance
                        break
                    end
                end
                %ev is produced after contraction steps
                %tolm will also be used below
            
                % SECTION 2: NEWTON-KANTOROVICH ITERATIONS
                if opt.printfxp>0
                    fprintf('\n');
                    fprintf('Begin Newton-Kantorovich iterations (for the %d. time)\n',k);
                    fprintf('  nwt          tol   \n');
                    
                end

                %do initial contraction iteration which is part of first N-K iteration
                [ev1, p1] = EntryExit.bellman(ev, Pi, P, beta,bparam); %also return choice probs=function of ev
                
                for nwt=1:opt.nstep; %do at most nstep of N-K steps

                    NKIter=NKIter+1;

                    %Do N-K step
                    dev = EntryExit.dbellman(p1, P, beta);
                    F = speye(n_state*2) - dev;
                    ev=ev-F\(ev-ev1); %resuing ev here
                    
                    %Contraction step for the next N-K iteration
                    [ev2, p1]= EntryExit.bellman(ev, Pi, P, beta,bparam); %also return choice probs=function of ev

                    %Measure the tolerance
                    tolm1=max(max(abs(ev-ev2)));

                    if opt.printfxp>0
                        fprintf('   %d     %16.8e  \n',nwt,tolm1);
                    end
                    if opt.printfxp>1
                        %plot ev
                        hold on
                        subplot(2,1,2), plot(ev -ev(1), '-r');
                    end

                    %Discard the N-K step if tolm1 got worse AND we can still switch back to SA
                    if (tolm1/tolm > 1+opt.rtolnk) && (k<opt.max_fxpiter);
                        if opt.printfxp>0
                            %Discrading the N-K step
                            fprintf('Discarding N-K step\n');
                        end
                        ev=ev1; %new contraction step should start from ev1
                        break;
                    else
                        ev1=ev2; %accept N-K step and prepare for new iteration
                    end;

                    %adjusting the N-K tolerance to the magnitude of ev
                    adj=ceil(log(abs(max(max(ev1)))));
                    ltol=opt.ltol0*10^adj;  % Adjust final tolerance
                    ltol=opt.ltol0;

                    if (tolm1 < ltol);
                        %N-K converged 
                        converged=true;
                        break
                    end

                end %Next N-K iteration

                if converged
                    if opt.printfxp>0
                        fprintf('Convergence achieved!\n\n');
                    end
                    break; %out of poly-algorithm loop
                else
                    if nwt>=opt.nstep
                        warning('No convergence! Maximum number of iterations exceeded without convergence!');
                        break; %out of poly-algorithm loop with no convergence
                    end
                end
            end
        end  % End of EntryExit.solve
        
        
        function theta_hat =  SequentialEstimationEE(datasim,theta_vec0,S)
            opt.tol  = 1e-9;
            opt.max_iter = 100;
            opt.output = 0;
            diff      = 1;
            iter      = 0;
            
            p0 = 0.5 * ones(S.n_state,1);
            ev0 = zeros(S.n_state,1);
            while(diff > opt.tol)
                f = @(theta_vec) - ImportDecision.ll_NFXP(theta_vec,datasim,ev0,S);
                Prob = conAssign(f,[],[],[],[],[],'2typeMixtureDDC',theta_vec0);
                Result = tomRun('Snopt',Prob);
                theta_vec = Result.x_k;
                
                if opt.output > 0
                    fprintf('Iteration: %d, Difference: %.4f, Time elapsed: %f Seconds \n',iter,diff,Result.REALtime);
                end                %STEP3: Update the weight using the theta
                
                [ll,p1] = ImportDecision.ll_NFXP(theta_vec,datasim,ev0,S);
                diff = max(abs(p1 - p0));
                iter = iter + 1;
                p0 = p1;
                theta_vec0 = theta_vec;
                if iter > opt.max_iter 
                    break
                end
                
            end
        end %End of Sequential EE
        
        
        
    end % End of static methods
end % End of class
