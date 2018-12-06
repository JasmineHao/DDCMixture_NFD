classdef EulerEquation
%FUNCTION_NAME - One line description of what the function or script performs (H1 line)
%Optional file header info (to give more details about the function than in the H1 line)
%Optional file header info (to give more details about the function than in the H1 line)
%Optional file header info (to give more details about the function than in the H1 line)
%
% Syntax:  [output1,output2] = function_name(input1,input2,input3)
%
% Inputs:
%    input1 - Description
%    input2 - Description
%    input3 - Description
%
% Outputs:
%    output1 - Description
%    output2 - Description
%
% Example: 
%    Line 1 of example
%    Line 2 of example
%    Line 3 of example
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: OTHER_FUNCTION_NAME1,  OTHER_FUNCTION_NAME2

% Author: FirstName FamilyName
% Work address
% email: 
% Website: http://www.
% May 2004; Last revision: 12-May-2004

%------------- BEGIN CODE --------------
    properties
    datasim; %the data (a,x,y);
    state;
    S;
    end

    methods (Static)
        

        
       %% ----------------------------------------------------
        %  The ee likelihood
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

        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%The past is in the past
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [logl, p1] =  likelihood(theta_vec,datasim,p1,S,InvF0,F_til)
            % GMM to update theta using Euler Equation
            % Used in two-step estimator of theta
            % This function gives sum of log-likelihood when fixing the p1(h-1)
            a = reshape(datasim.at,S.N*S.T,1);
            z = reshape(datasim.zt,S.N*S.T,1);
            y = reshape(datasim.yt,S.N*S.T,1);
            [n_state,d_state] = size(S.state);
            ind = (y - 1) * n_state + z ;
            ccp = reshape(p1,n_state,2);
            %
            dvpdth   = exp(S.state(:,5)).* [ones(n_state,1),S.state(:,1),S.state(:,2)]; %Derivative of variable profit
            dfcdth   = [ones(n_state,1),S.state(:,3)];% Derivative of fixed cost
            decdth   = [ones(n_state,1),S.state(:,4)];
            % Step 2: Compute the derivative of contraction mapping w.r.t mp
            dpidth_1 =  [dvpdth,-dfcdth,zeros(size(decdth))];
            dpidth_0 =  [dvpdth,-dfcdth,-decdth];
            pi       = [dpidth_1 * theta_vec ;dpidth_0 * theta_vec] ;
            p0       = 1 - p1; %The p_0 from last period
            V        = InvF0 * ( 0.5772156649 - log(p0) ) ;
          %
            tmp = pi + S.beta * F_til * V;
            p1_1 = exp(tmp) ./ ( 1 + exp(tmp));
            data_p1 = p1_1(ind);
            %
            logl  = sum([a==1].*log(data_p1) + [a==2].*log(1 - data_p1));
            if nargout > 1
                p1 = reshape(p1_1,n_state * 2,1);
            end
        end

        function g_logl =  g_likelihood(theta_vec,datasim,p1,S)
            % GMM to update theta using Euler Equation
            % Provide the gradient of loglikelihood function
            a = reshape(datasim.at,S.N*S.T,1);
            z = reshape(datasim.zt,S.N*S.T,1);
            y = reshape(datasim.yt,S.N*S.T,1);
            [n_state,d_state] = size(S.state);
            ind = (y - 1) * n_state + z ;
            % Reshape the estimation of ccp from last iteration
            ccp  = reshape(p1,n_state,2);
            dvpdth   = exp(S.state(:,5)).* [ones(n_state,1),S.state(:,1),S.state(:,2)]; %Derivative of variable profit
            dfcdth   = [ones(n_state,1),S.state(:,3)];% Derivative of fixed cost
            decdth   = [ones(n_state,1),S.state(:,4)];
            % Step 2: Compute the derivative of contraction mapping w.r.t mp
            dpidth_1 =  [dvpdth,-dfcdth,zeros(size(decdth))];
            dpidth_0 =  [dvpdth,-dfcdth,-decdth];
            dpidth   = [dpidth_1;dpidth_0];
            pi       = dpidth * theta_vec;
            %
            v_til_x  = pi + repmat(S.beta * S.P *( log(1- ccp(:,2)) - log(1 - ccp(:,1)) ),2,1);
            p1_1     = exp(v_til_x) ./ ( 1 + exp(v_til_x));
            p1_1     = reshape(p1_1,n_state,2);
            data_p1 = p1_1(ind);
            g_logl  = sum(bsxfun(@times,([a==1] - data_p1),dpidth(ind,:)));
        end
        
        
        function[v_til,p1] = solve(P,state,beta,theta_vec)
            %----------------------------------------------------------------------
            %  SYNTAX: [ev,p1] = EulerEquation.solve(P,state,beta,theta);
            %----------------------------------------------------------------------
            %  INPUT:
            %     EV0 :      m x 1 matrix or scalar zero. Initial guess of choice specific expected value function, EV.
            %     param:     parameters
            %     state:     State transition matrix
            %     P:         State transition of exogeneous variable
            %----------------------------------------------------------------------
            %
            %  OUTPUT:
            %     P1:     m x 1 matrix of conditional choice probabilities, P(d=enter|x)
            %     V_TIL:  m x 1 matrix of difference in value functions
            %             v(1,x_t) - v(0,x_t)
            %     F:      m x m matrix of derivatives of Identity matrix I minus
            %             contraction mapping operator, I-T' where T' refers to derivative of the expected  value function
            %----------------------------------------------------------------------
            opt.max_cstp   = 30;
            opt.min_cstp   = 4;
            opt.printfxp   = 0;
            opt.ctol       = 1e-8;
            %opt.method     = 'value';
            opt.method     = 1;
            n_state = size(state,1);

            dvpdth   = exp(state(:,5)).* [ones(n_state,1),state(:,1),state(:,2)]; %Derivative of variable profit
            dfcdth   = [ones(n_state,1),state(:,3)];% Derivative of fixed cost
            decdth   = [ones(n_state,1),state(:,4)];
            % Step 2: Compute the derivative of contraction mapping w.r.t mp
            dpidth_1 =  [dvpdth,-dfcdth,zeros(size(decdth))];
            dpidth_0 =  [dvpdth,-dfcdth,-decdth];
            pi       = [dpidth_1 * theta_vec ,dpidth_0 * theta_vec] ;

             %Difference in one period payoff
            tolm = 1;
            if opt.method > 0
                v_til = zeros(n_state,2); %Initial guess of v_til
                for j=1:opt.max_cstp
                    v_til1 = EulerEquation.contraction_v(v_til, pi , beta, P);
                    tolm1=max(max(abs(v_til1 - v_til)));
                    rtolm=tolm1/tolm;
                    if opt.printfxp>0
                         fprintf(' %3.0f   %16.8f %16.8f\n',j, tolm1,rtolm);
                    end
                                            %prepare for next iteration
                    v_til=v_til1;
                    tolm=tolm1;

                        %stopping criteria
                    if (j>=opt.min_cstp) && (tolm1 < opt.ctol)
                            %go to NK iterations due to absolute tolerance
                         break;
                    end
                end
                p1 = 1 ./ (exp(-v_til) + 1);
                p1    = reshape(p1,n_state,2);
                v_til = reshape(v_til,n_state,2);
            else
                p_10  = 0.5 * ones(n_state,2);
                for j = 1:opt.max_cstp
                    p_1 = EulerEquation.contraction_p(p_10,pi,beta,P);
                    tolm1=max(max(abs(p_1 - p_10)));
                    rtolm=tolm1/tolm;
                    if opt.printfxp>0
                         fprintf(' %3.0f   %16.8f %16.8f\n',j, tolm1,rtolm);
                    end
                    p_10 = p_1;
                    tolm=tolm1;
                        %stopping criteria
                    if (j>=opt.min_cstp) && (tolm1 < opt.ctol)
                            %go to NK iterations due to absolute tolerance
                         break;
                    end
                end
                p1    = reshape(p_1,n_state,2);
                v_til = log(p_1 ./ (1 - p_1));
            end
        end


        function v_til = contraction_v(v_til, pi , beta, P)
        %------------------------------------------------------------------
        % contraction mapping of Euler Equation
        %------------------------------------------------------------------
        n_state = size(P,2);
        pi   = reshape(pi,n_state,2);
        A = 0.5 * beta * P * (1 ./ (1 + exp(v_til(:,1))) + ...
            1 ./ (1 + exp(v_til(:,1)) ) ) .* (pi(:,1) - pi(:,2) + v_til(:,2) - v_til(:,1) );
        B = beta * P * (log(1+exp(v_til(:,1))) - log(1+exp(v_til(:,2))) );
        v_til(:,1) = pi(:,1) + B ;
        v_til(:,2) = pi(:,2) + B  ;
        end

        function p1 = contraction_p(p10,pi,beta,P)
            p10 = 1 - p10;
            tmp1 = pi(:,1) + beta * P *( log(p10(:,2)) - log(p10(:,1)) );
            tmp2 = pi(:,2) + beta * P *( log(p10(:,2)) - log(p10(:,1)) );
            p1(:,1) = exp(tmp1) ./ ( 1 + exp(tmp1));
            p1(:,2) = exp(tmp2) ./ ( 1 + exp(tmp2));
        end

        function err = constraint(p1,theta_vec,S)
            n_state = size(S.P,2);

            dvpdth   = exp(S.state(:,5)).* [ones(n_state,1),S.state(:,1),S.state(:,2)]; %Derivative of variable profit
            dfcdth   = [ones(n_state,1),S.state(:,3)];% Derivative of fixed cost
            decdth   = [ones(n_state,1),S.state(:,4)];
            % Step 2: Compute the derivative of contraction mapping w.r.t mp
            dpidth_1 =  [dvpdth,-dfcdth,zeros(size(decdth))];
            dpidth_0 =  [dvpdth,-dfcdth,-decdth];
            pi       = [dpidth_1 * theta_vec ,dpidth_0 * theta_vec] ;

            ccp =  reshape(p1,n_state,2);
            p1  = reshape(p1,n_state*2,1);
            tmp1 = pi(:,1) + S.beta * S.P *( log(1 - ccp(:,2)) - log(1 - ccp(:,1)) );
            tmp2 = pi(:,2) + S.beta * S.P *( log(1 - ccp(:,2)) - log(1 - ccp(:,1)) );
            tmp = [tmp1;tmp2];
            err = (exp(tmp)./(1+exp(tmp)) - p1);
        end
        
        function [theta_vec,p1] = estimation_2step(theta_vec0,datasim,S)
            opt.max_cstp   = 30;
            opt.min_cstp   = 2;
            opt.printfxp   = 1;
            opt.ctol       = 1e-7;
            Name = 'The Problem';

            p10 = 0.5 * ones(size(S.P,2)*2,1);
            tolm = 1;
            F_0 = [zeros(size(S.P)), S.P;zeros(size(S.P)), S.P];
            InvF0 = inv(eye(S.n_state * S.n_action) - S.beta * F_0);
            F_til = [S.P,-S.P;S.P,-S.P];
            for i = 1: opt.max_cstp
                f = @(theta_vec)(-EulerEquation.likelihood(theta_vec,datasim,p10,S,InvF0,F_til));
                g = @(theta_vec)(-EulerEquation.g_likelihood(theta_vec,datasim,p10,S));
                Prob = conAssign(f,g ,[],[],[],[],Name,theta_vec0);
                Result = tomRun('Snopt',Prob);
                theta_vec = Result.x_k;

                [logl,p1] = EulerEquation.likelihood(theta_vec,datasim,p10,S,InvF0,F_til);
                tolm1 = max(max(abs(p10 - p1)));
                rtolm = tolm1 / tolm;
                p10 = p1;
                theta_vec0 = theta_vec;
                if opt.printfxp > 0
                    fprintf(' %3.0f   %16.8f %16.8f %16.8f\n',i, tolm1,rtolm,logl);
                end
                tolm = tolm1;
                if (tolm < opt.ctol ) && (i > opt.min_cstp)
                    break
                end
            end
        end

        function [theta_vec,p1] = estimation_mpec(theta_vec0,datasim,p1,S)
            n_state = size(S.state,1);
            theta_p_vec = [theta_vec0;p1] ;
            cons = @(theta_p_vec) EulerEquation.constraint(theta_p_vec(8:end),theta_p_vec(1:7), S);
            obj = @(theta_p_vec)(-EulerEquation.likelihood(theta_p_vec(1:7),datasim,theta_p_vec(8:end),S));
            g = @(theta_p_vec)(-EulerEquation.g_likelihood(theta_p_vec(1:7),datasim,theta_p_vec(8:end),S));
            Name = 'The Problem';
            x_L = [-inf * ones(7,1);zeros( 2 * n_state,1 )];
            x_U = [inf * ones(7,1);ones( 2 * n_state,1 )];
            c_L = 1e-2 * ones(2*n_state,1);
            c_U = 1e-2 * ones(2*n_state,1);

            Prob = conAssign(obj,g,[],[],x_L,x_U ,Name,theta_p_vec, ...
                              [] , [], [],[], [], cons, [], [], [], c_L, c_U);
            Result = tomRun('Knitro',Prob);
            theta_vec = Result.x_k(1:7);
            p1  = Result.x_k(8:end);
            PrintResult(Result,2);;
        end %End of Estimation using MPEC
    end %End of Methods
    
end
