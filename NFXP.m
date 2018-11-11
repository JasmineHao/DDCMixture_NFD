classdef NFXP
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

        %------------------------------------------------------------------
        function pi = profit(theta_vec,S)
            n_state = S.n_state;
            dvpdth   = exp(S.state(:,5)).* [ones(n_state,1),S.state(:,1),S.state(:,2)]; %Derivative of variable profit
            dfcdth   = [ones(n_state,1),S.state(:,3)];% Derivative of fixed cost
            decdth   = [ones(n_state,1),S.state(:,4)];
            % Step 2: Compute the derivative of contraction mapping w.r.t mp
            dpidth_1 =  [dvpdth,-dfcdth,zeros(size(decdth))];
            dpidth_0 =  [dvpdth,-dfcdth,-decdth];
            pi       = [dpidth_1 * theta_vec , dpidth_0 * theta_vec] ;
        end %end of profit

        %%%----------------------------------------------------------------
        function [ev1, p1]=bellman(ev, pi, S)
            n_action = 2;
            n_state = S.n_state;
            ev = reshape(ev,n_state, n_action);
            ev_11 = pi(:,1) + S.beta * S.P * ev(:,1);
            ev_01 = pi(:,2) + S.beta * S.P * ev(:,1);
            ev_00  = S.beta * S.P * ev(:,2);
            ev_1= max(ev_11,ev_00);
            ev_0= max(ev_01,ev_00);
            ev1 = [ev_1;ev_0];
            if nargout>1
                p1 = [1 ./ ( 1 + exp(ev_00 - ev_11)); 1 ./ ( 1 + exp(ev_00 - ev_01))];
            end
        end

        %% %---------------------------------------------------------------

        function dev = dbellman(p1, S)
            % Bellman equation
            % NFXP.DBELMANN:   Procedure to compute Frechet derivative of Bellman operator
            % Syntax :          ev=nfpx.dbellman(pk,  P, mp);
            n_state = S.n_state;
            p1 = reshape(p1,n_state,2);
            %tmp11 : Prob(a=1|y=1) , tmp10 : Prob(a=0|y=1)
            tmp11= S.P .* repmat(p1(:,1),1,n_state); %y = 1, a = 1
            tmp10= S.P .* repmat((1 - p1(:,1)),1,n_state); %y = 1, a = 0
            tmp01= S.P .* repmat(p1(:,2),1,n_state);%y = 0, a = 1
            tmp00= S.P .* repmat((1 - p1(:,2)),1,n_state);%y = 0, a = 0
            dev   = S.beta * [tmp11 , tmp10 ; tmp01, tmp00];
        end

        %%-----------------------------------------------------------------
        function [ev,p1,F] = solve(theta_vec,S)
            %--------------------------------------------------------------
            %SYNTAX: [ev,p1,F] = NFXP.solve(theta_vec,S)
            %--------------------------------------------------------------
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
            n_action = S.n_action;
            n_state  = S.n_state;
            pi       = NFXP.profit(theta_vec,S);

            % Initial guess of ev
            ev = zeros(n_state * n_action, 1 );
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
                for j=1:opt.max_cstp
                    ev1=NFXP.bellman(ev, pi, S);

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
                    if (j>=opt.min_cstp) && (abs(S.beta-rtolm) < opt.rtol)
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
                [ev1, p1] = NFXP.bellman(ev, pi, S); %also return choice probs=function of ev

                for nwt=1:opt.nstep %do at most nstep of N-K steps

                    NKIter=NKIter+1;

                    %Do N-K step
                    dev = NFXP.dbellman(p1, S);
                    F = speye(n_state*2) - dev;
                    ev=ev-F\(ev-ev1); %resuing ev here

                    %Contraction step for the next N-K iteration
                    [ev2, p1]= NFXP.bellman(ev, pi, S); %also return choice probs=function of ev

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
        end  % End of NFXP.solve
 
        function err = constraint(theta_vec,ev,S)
            %This one use the same notation as EulerEquation
            n_state = size(S.P,2);
            ev_0 = reshape(ev,n_state,2);
            pi       = NFXP.profit(theta_vec,S) ;
            tmp1     = [pi(:,1) + S.beta * S.P * ev_0(:,1); pi(:,2) + S.beta * S.P * ev_0(:,1)];
            tmp2     = repmat(S.beta * S.P * ev_0(:,2),2,1);
            err      = max(tmp1,tmp2) - ev;
        end

        function d_cons = d_constraint(theta_vec,ev,S)
            n_state = size(S.P,2);
            ev_0 = reshape(ev,n_state,2);
            dvpdth   = exp(S.state(:,5)).* [ones(n_state,1),S.state(:,1),S.state(:,2)]; %Derivative of variable profit
            dfcdth   = [ones(n_state,1),S.state(:,3)];% Derivative of fixed cost
            decdth   = [ones(n_state,1),S.state(:,4)];
            % Step 2: Compute the derivative of contraction mapping w.r.t mp
            dpidth_1 =  [dvpdth,-dfcdth,zeros(size(decdth))];
            dpidth_0 =  [dvpdth,-dfcdth,-decdth];
            dpidth   =  [dpidth_1;dpidth_0];
            pi       = [dpidth_1 * theta_vec , dpidth_0 * theta_vec] ;
            tmp1     = [pi(:,1) + S.beta * S.P * ev_0(:,1); pi(:,2) + S.beta * S.P * ev_0(:,1)];
            tmp2     = repmat(S.beta * S.P * ev_0(:,2),2,1);
            d  = tmp1 >= tmp2;
            dcondth = diag(d) * dpidth;
            d  = reshape(d,n_state,2);
            dcondev = S.beta * [diag(d(:,1)) * S.P , diag(1 - d(:,1)) * S.P ;
                diag(d(:,2)) * S.P,diag(1-d(:,2)) * S.P  ] - eye(n_state * 2);
            d_cons = [dcondth,dcondev];
        end
        function logl = likelihood(theta_vec,ev,datasim,S)
            %-------------------------------------------------------------
            %--Likelihood function
            %-------------------------------------------------------------
            a = reshape(datasim.at,S.N*S.T,1);
            z = reshape(datasim.zt,S.N*S.T,1);
            y = reshape(datasim.yt,S.N*S.T,1);
            [n_state,d_state] = size(S.state);
            ind = (y - 1) * n_state + z ;

            dvpdth   = exp(S.state(:,5)).* [ones(n_state,1),S.state(:,1),S.state(:,2)]; %Derivative of variable profit
            dfcdth   = [ones(n_state,1),S.state(:,3)];% Derivative of fixed cost
            decdth   = [ones(n_state,1),S.state(:,4)];
            % Step 2: Compute the derivative of contraction mapping w.r.t mp
            dpidth_1 =  [dvpdth,-dfcdth,zeros(size(decdth))];
            dpidth_0 =  [dvpdth,-dfcdth,-decdth];
            pi       = [dpidth_1 * theta_vec ,dpidth_0 * theta_vec] ;

            ev = reshape(ev,n_state,2);

            ev_11 = pi(:,1) + S.beta * S.P * ev(:,1);
            ev_01 = pi(:,2) + S.beta * S.P * ev(:,1);
            ev_00  = S.beta * S.P * ev(:,2);
            p1 = [1 ./ ( 1 + exp(ev_00 - ev_11)); 1 ./ ( 1 + exp(ev_00 - ev_01))];
            data_p1 = p1(ind);
            logl  = sum([a==1].*log(data_p1) + [a==2].*log(1 - data_p1));
        end
            % Compute score
        function g = g_likelihood(theta_vec,ev,datasim,S)
            a = reshape(datasim.at,S.N*S.T,1);
            z = reshape(datasim.zt,S.N*S.T,1);
            y = reshape(datasim.yt,S.N*S.T,1);
            [n_state,d_state] = size(S.state);
            ind = (y - 1) * n_state + z ;
            ev_0 = reshape(ev,n_state,2);
            dvpdth   = exp(S.state(:,5)).* [ones(n_state,1),S.state(:,1),S.state(:,2)]; %Derivative of variable profit
            dfcdth   = [ones(n_state,1),S.state(:,3)];% Derivative of fixed cost
            decdth   = [ones(n_state,1),S.state(:,4)];
            % Step 2: Compute the derivative of contraction mapping w.r.t mp
            dpidth_1 =  [dvpdth,-dfcdth,zeros(size(decdth))];
            dpidth_0 =  [dvpdth,-dfcdth,-decdth];
            dpidth   = [dpidth_1;dpidth_0];

            pi       = [dpidth_1 * theta_vec ;dpidth_0 * theta_vec] ;
            v1       = pi + S.beta * repmat(S.P * ev_0(:,1),2, 1);
            p1      = exp(v1) ./ ( 1 + exp(v1));
            p1 = reshape(p1,n_state,2);

            dtdth_1  =  p1(:,1) .* dpidth_1 ; %dpi/dtheta when y = 1
            dtdth_0  =  p1(:,2) .* dpidth_0 ; %dpi/dtheta when y = 0
            dtdth    = [dtdth_1;dtdth_0];
            dpidth   = [dpidth_1;dpidth_0];
                %tmp11 : Prob(a=1|y=1) , tmp10 : Prob(a=0|y=1)
                tmp11= S.P .* repmat(p1(:,1),1,n_state); %y = 1, a = 1
                tmp10= S.P .* repmat((1 - p1(:,1)),1,n_state); %y = 1, a = 0
                tmp01= S.P .* repmat(p1(:,2),1,n_state);%y = 0, a = 1
                tmp00= S.P .* repmat((1 - p1(:,2)),1,n_state);%y = 0, a = 0
                dev   = S.beta * [tmp11 , tmp10 ; tmp01, tmp00];
                F = speye(n_state*2) - dev;
                dvdth   = F \ dtdth;
                dvdth_1 = dvdth(1:n_state,:);
                dvdth_0 = dvdth((1+n_state):(2*n_state),:);
                devdth  = S.P * ( dvdth_1 - dvdth_0 );
                g =  sum(bsxfun(@times,( (a == 1) - p1(ind)) , dpidth(ind,:) +...
                  S.beta *  devdth(z,:)));
        end




        function [thetahat, thetase, thetacov, g, ll, iterinfo]=estimation_newton(llfun, theta_vec0);
            %%
            %----------------------------------------------------------------------
            %   NFXP.MAXLIK:  Routine to maximize likelihood functions based on
            %               - Analytical gradients.
            %               - BHHH Algorithm
            %               - Line search (step-halfing)
            %   SYNTAX:  = NFXP(datasim,theta,S)
            %----------------------------------------------------------------------
            %
            %   OUTPUTS:
            %       thetahat   :   Estimated parameters (structure with as many fields as p0)
            %       thetase    :   Standard errors of estimated parameters (structure with number of fields equal to number elements in pnames, nP)
            %       thetacov   :   Variance Covariance Matrix (nP x nP matrix)
            %       g          :   Score vectors (N x nP)
            %       ll         :   Likelihood contributions (N x 1 vector)
            %       convflag   :   Bolean (true/false) idicator of convergence
            %----------------------------------------------------------------------
            %   INPUTS:
            %       llfun      :  Function handle for likelihood function
            %
            %       theta0     :  Structure holding starting values of parameters
            %
            %       pnames     :  k dimensional cell array, that holds names of the field in the parameter structure
            %                     to be estimated. Then likelihood is only maximized with respect to the parameters in pnames
            %                     and it is thus possible to fix parameters at the values set in theta0
            %       options    :  Structure holding voluntary options for ML. If unspecified default options are used
            %                     options structure has the following fields:
            %
            %           options.title (String):
            %                   A title for the model.
            %
            %           options.output (Scalar) :
            %                   0 if no output, 1 for final output, 2 for detailed iteration output (default is 2)
            %
            %           options.maxit (Scalar):
            %                   Maximum number of iterations (default is 100)
            %
            %           options.tolgdirec (Scalar):
            %                   Tolerance for gradient'direc (default is 1e-3)
            %
            %           options.lambda0:
            %                   The basic step length (default is 1).
            %----------------------------------------------------------------------

            tt=tic;
            par.title={''};
            par.maxit=100;
            par.tolgdirec=1e-2;
            par.tol_bhhh=0.1;
            par.output=2;
            par.lambda0=1;
            par.hess=0;
            theta_vec1=theta_vec0;

            global BellmanIter NKIter
            iterinfo.MajorIter=0;
            iterinfo.ll=0;
            iterinfo.NKIter=0;
            iterinfo.BellmanIter=0;

            if nargin==4
                pfields=fieldnames(options);
                for i=1:numel(pfields)
                    par.(pfields{i})=options.(pfields{i});
                end
            end

            lambda=par.lambda0;        % initialize stepsize
            k     = length(theta_vec0);     % Number of parameters
            % -------------------------
            %  ** BHHH Optimization **
            % -------------------------

%             if (par.hess==1)
%                 strhess='Analytical Hessian';
%                 [ll0,ev,s0, h]=llfun(theta0);
%                 h=-h;
%             elseif (par.hess==0);
%                 strhess='Outer product of scores';
%                 [ll0,ev,s0]=llfun(theta0);
%                 h=s0'*s0;
%             end;
%%
            strhess='Outer product of scores';
            [ll0,ev,s0]=llfun(theta_vec0);
            h=s0'*s0;

            iterinfo.ll=iterinfo.ll+1;
            iterinfo.NKIter=iterinfo.NKIter+NKIter;
            iterinfo.BellmanIter=iterinfo.BellmanIter+BellmanIter;


            g=sum(s0)';
            direc=h\g;
            tolm=g'*direc;

            %%
            if tolm<0 % convex area of likelihood, Newton-Raphson moves away from max
                h=s0'*s0;
                direc=h\g;
                tolm=g'*direc;
                if par.output>=1
                    disp('Convex area of likelihood: Switch to BHHH')
                end
            end
            %%
            iterinfo.MajorIter=0;
            for it=0:par.maxit
                for i=1:k
                    theta_vec1(i) =theta_vec0(i) + lambda*direc(i);
                end
                ev1 = ev;
                if lambda==par.lambda0
                % If previous step was accepted
%                      if ((par.hess==0) | (abs(tolm) > par.tol_bhhh));
%                         strhess='Outer product of the scores';
%                         [ll,ev,s]=llfun(p1);
%                         h=s'*s;
%                     else
%                         [ll,ev,s, h]=llfun(p1);
%                         strhess='Analytical Hessian';
%                         h=-h;
%                     end;
                    [ll,ev,s]=llfun(theta_vec1);
                    h=s'*s;
                    iterinfo.ll=iterinfo.ll+1;
                    iterinfo.NKIter=iterinfo.NKIter+NKIter;
                    iterinfo.BellmanIter=iterinfo.BellmanIter+BellmanIter;
                else
                    [ll, ev, s]=llfun(theta_vec1); % Do we need to recalculate s, when line searching
                    iterinfo.ll=iterinfo.ll+1;
                    iterinfo.NKIter=iterinfo.NKIter+NKIter;
                    iterinfo.BellmanIter=iterinfo.BellmanIter+BellmanIter;
                end

                if sum(ll)<sum(ll0)
                    lambda=lambda/2;
                    if par.output>1
                        fprintf('\n');
                        fprintf('.............LINESEARCHING %f %f \n',sum(ll0),sum(ll));
                        fprintf('\n');
                    end
                    if lambda<=.01
                        if par.output>1
                           % fprintf('.............WARNING: Failed to increase - downhill step accepted\n');
                        end
%                         theta0=theta1;
%                         ll0=ll;                 % Accept step
%                         lambda=par.lambda0;
                    end
                else       % IF INCREASE IN LIKELIHOOD, CALCULATE NEW DIRECTION VECTOR
                    theta_vec0=theta_vec1;
                    ll0=ll;                 % Accept step


                    % ** Plot iteration info ***
                    if par.output>1
                        fprintf('Iteration: %d   \n', iterinfo.MajorIter);
                        fprintf('grad*direc       %10.5f \n', g'*direc);
                        fprintf('||rel grad||     %10.5f \n', max(abs(g/sum(ll0))));
                        fprintf('Log-likelihood   %10.5f \n', sum(ll));
                        fprintf('Hessian update   %10s \n', strhess);
                        fprintf('Step length      %10.5f \n', lambda);
                        fprintf('---------------------------------------------------------------------\n');
                        fprintf('        Estimates      Direction         Gradient\n');
                        i=0;
                        for i=1:k
                            fprintf('     %10.4f    %10.4f        %10.4f\n', theta_vec0(i), direc(i), g(i));
                        end
                        fprintf('--------------------------------------------------------------------\n');
                        fprintf('\n');
                    end
                    g=sum(s)';
                    iterinfo.MajorIter=iterinfo.MajorIter+1;
                    direc=h\g;
                    tolm=g'*direc;
                    tolgrad=max(abs(g/sum(ll0)));

                    if tolm<0 % convex area of likelihood, Newton-Raphson moves away from max
                        h=s'*s;
                        direc=h\g;
                        if par.output>=0
                            disp('Convex area of likelihood: Switch to BHHH')
                        end
                    end
                    lambda=par.lambda0; % and reset stepsize
                end

                if tolm < par.tolgdirec % Stopping rule
                    break;
                end

            end

            thetacov=inv(h);
            se=(sqrt(diag(thetacov)));
            thetahat=theta_vec0;
            for i=1:k
                thetase(i)=se(i);
            end

            % ---------------------------
            %  ** Plot final output  ***
            % ---------------------------
            if it<par.maxit
                if par.output >= 1
                    cr=corrcov(thetacov);
                    fprintf('---------------------------------------------------------------------\n');
                    fprintf('***                   Convergence Achieved                        ***\n');
                    fprintf('---------------------------------------------------------------------\n');
                    disp('                         _ ');
                    disp('                         \`\ ');
                    disp('                         |= | ');
                    disp('                        /-  ;.---. ');
                    disp('                  _ __.''     (____) ');
                    disp('                   `         (_____) ');
                    disp('                   _''  ._ .'' (____) ');
                    disp('                    `        (___) ');
                    disp('                   --`''------''` ');
                    fprintf('%s \n', char(par.title));
                    fprintf('Number of iterations: %d\n', iterinfo.MajorIter);
                    fprintf('grad*direc       %10.5f \n', g'*direc);
                    fprintf('Log-likelihood   %10.5f \n', sum(ll));
                    fprintf('Step length      %10.5f \n', lambda);
                    fprintf('\n');
                    fprintf('    %13s %13s %13s\n','Estimates','s.e.','t-stat');

                    fprintf('---------------------------------------------------------------------\n');
                    for i=1:k
                        fprintf('     %13.4f %13.4f %13.4f\n',  theta_vec0(i), se(i), theta_vec0(i)/se(i));
                        thetase(i) =se(i);
                    end
                    fprintf('---------------------------------------------------------------------\n');

                    fprintf('\n');

                    fprintf('Correlation matrix of parameters\n');
                    fprintf('%18s','');

                    if ismac & par.output >= 3
                        !say 'Success! Convergence Achieved'
                    end
                    tt_end=toc(tt);
                    fprintf('\n');
                    fprintf('Time to convergence is %3.0f min and %2.2f seconds\n', floor((tt_end)/60),tt_end - 60*floor((tt_end)/60));

                    fprintf('\n');

                end
                % output a flag that we did converge
                iterinfo.convflag = true;
            else
                fprintf('----------------------------------------------------------------------------\n');
                fprintf('***   BHHH failed to converge: Maximum number of iterations is reached *** \n');
                fprintf('----------------------------------------------------------------------------\n');
                disp('                   _,....._ ')
                disp('                  (___     `''-.__ ')
                disp('                 (____ ')
                disp('                 (____ ')
                disp('                 (____         ___ ')
                disp('                      `)   .-''` ')
                disp('                      /  .'' ')
                disp('                     | =| ')
                disp('                      \_\ ')

                if ismac
                    !say 'B triple H failed to converge: Maximum number of iterations is reached without convergence. Better luck next time'
                end

                % no luck this time
                iterinfo.convflag = false;

            end
        end % end of estimation_newton

    
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
            
            [devdth,dpidth] = NFXP.devdth(S,F,p1);
            % Step 4: Compute the derivative of log-likelihood w.r.t theta
            score    = bsxfun(@times,( (a == 1) - p1(ind)) , dpidth(ind,:) +...
                      S.beta *  devdth(z,:));  
        end
    
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

    end
end
