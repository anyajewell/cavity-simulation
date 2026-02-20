classdef OptClass < handle
    %OPTCLASS 
    
    properties
        m %number of inputs
        p %number of outputs
        mugap %gap in identification model F coefficients
        pt %number of tracking outputs
        nhat %identification order
        ell %MPC prediction horizon 
        
        Rkl %Reference command for tracking measurements (C_t) over horizon: num tracking measurements x mpc horizon
        %in the paper it is num tracking measurements*mpc horizon x 1, but
        %we use the matrix format to make updating a little easier
        
        Fbar %matrix of model denominator coefficents ( p*ell x p*(n+ell) )
        Gbar %matrix of model numerator coefficents ( p*ell x m*(n+ell) ) 
        Gbarf %matrix of model numerator coefficents ( p*ell x 1*(n+ell) )
        Gbarc %matrix of model numerator coefficents ( p*ell x 1*(n+ell) )
        
        %Observability-like and Toeplitz matrices for response prediction
        Gamma %Matrix made from past y and u data ( p*ell x 1 )
        T     %Toeplitz matrix (p*ell x m*ell)
        Tfi    %Toeplitz matrix (p*ell x m*ell)
        Tco    %Toeplitz matrix (p*ell x m*ell)
        
        %Constraint parameters
        %C_t,C_c are matrices that "pick" which measurements we want
        %tracked or included in the constraints at a single step (y_t,k = C_t*y,k). 
        %The following matrices do the same but for multiple measurements across the mpc horizon
        C_t
        C_tl %C_t but expanded for every measurement across horizon (mpc horizon * num tracked measurements x mpc horizon * num measurements)
        C_l %C_c but expanded for every measurement across horizon (mpc horizon * num measurement constraints x mpc horizon * num measurements)
        Dconstr %script D constraint from C_l*Y1|k <= Dconstr
        Delu_max %max change in control input, deltau
        Delu_min %min change in control input, deltau
        u_max %max control input, u
        u_min %min control input, u
        Aconstr %Constraint matrix for contraints of type Ax <= b
        bconstr %Constraint vector for constraints of type Ax <= b
         

        %Optimization Weights
        Q
        Qi
        Ru
        Rdu 
        S

        %Tracking Error Integrator State
        i
        
        %deltaU creation matrices
        %TODO: Explain what these are with dimensions
        bc
        Shm 
        bcp
        ShpInv 
 
        
        %MPC computed controls 
        Umpc %the optimization variable across the horizon (num inputs x mpc horizon) 
        UC   %Coarse control sequence
        uOut %the control input we actually use (num inputs x 1) 
        slack %the resulting slack variable from the optimization
        
        %Optimization algorithm 
        optAlg %'quadprog' or 'grad' or 'admm'
        numIter %number of iterations (max iterations for quadprog, grad, and admm)
        lam %Lagrange multipliers for grad has dimensions of [Umpc;slack]
        
        %Option Flags
        slackOn %=1 if slack active - slack is ignored if admm is selected`
        numCon %Number of output constraints, = 0 if no output constraints 

        %equality constraint matrices
        beqDyn
        AeqDyn

        %flag for enabling explicit inverse
        EGT_flag 
    end
    
    methods (Access = public)
        function obj = OptClass(opts_, r_ , specialFlag_ ) 
            %OPTCLASS Construct an instance of this class 
            
            %Initialize sizes
            obj.m      = opts_.m; %number of inputs
            if specialFlag_
                obj.m = 2;
            end
            obj.p      = opts_.p; %number of outputs
            obj.nhat   = opts_.nhat; %Identification order
            obj.mugap  = opts_.mugap; %Identification F coefficient gap
            obj.ell    = opts_.ell; %MPC prediction horizon 
            
            %Initialize the Observability and Toeplitz matrices for
            %prediction  
            obj.Gamma = zeros(opts_.ell*opts_.p,1);


            obj.T     = zeros(opts_.ell*opts_.p,opts_.ell*opts_.m); 
            obj.Tfi    = zeros(opts_.ell*opts_.p,opts_.ell*1); 
            obj.Tco    = zeros(opts_.ell*opts_.p,opts_.ell*1); 
            
            
            %Initialize the Fbar Gbar matrices for GPC-like prediction
            obj.Fbar = sparse( zeros( opts_.p*opts_.ell , opts_.p*(opts_.nhat + opts_.ell ) ) );
            obj.Fbar( : , opts_.p*opts_.nhat + 1 : end ) = eye( opts_.p * opts_.ell );

            obj.Gbar  = sparse( zeros( opts_.p*opts_.ell , opts_.m*(opts_.nhat + opts_.ell ) ) ); 
            obj.Gbarf = sparse( zeros( opts_.p*opts_.ell , 1*(opts_.nhat + opts_.ell ) ) );
            obj.Gbarc = sparse( zeros( opts_.p*opts_.ell , 1*(opts_.nhat + opts_.ell ) ) );
            
            %Command to follow over horizon
            obj.Rkl = kron(ones(1,opts_.ell),r_); %r_*ones(size(opts_.C_t,1),opts_.ell);  
             
  
            %Initialize Tracking and Constraint matrices 
            obj.C_t  = opts_.C_t;
            obj.pt   = size(obj.C_t,1);
            obj.C_tl = kron(eye(opts_.ell), opts_.C_t );
            obj.C_l  = kron(eye(opts_.ell), opts_.C * opts_.C_c ); % C_temp is the script C from the constraint C*y_c + D < 0
            obj.Dconstr = opts_.Dconstr;
            obj.Delu_max = opts_.Delu_max;
            obj.Delu_min = opts_.Delu_min;
            obj.u_max = opts_.u_max;
            obj.u_min = opts_.u_min;
 
                 
            %Initialize optimization weights
            obj.Q    = blkdiag(opts_.Qhat, opts_.Phat);
            obj.Qi   = opts_.Qi*eye(obj.pt*opts_.ell);
            obj.Ru   = opts_.Ru;
            obj.Rdu  = opts_.Rdu; 
            obj.S    = opts_.S;

            

            %Tracking Error Integrator State
            obj.i    = zeros(obj.pt,1);
            











            %Set Option Flags
            if det( opts_.S ) == 0
                obj.slackOn = 0;
            else
                obj.slackOn = 1;
            end
            if strcmp(opts_.optAlg,'admm') 
                %since admm uses a special slack we disable PCAC slack
                obj.slackOn = 0;
            end

            % Number of constraints on yk
            obj.numCon     = size( opts_.C , 1 ); 

             
            %Initialize matrices for creating deltaU 
            obj.bc   = [-eye(opts_.m); zeros((opts_.ell-1) * opts_.m, opts_.m)]; 
            prototype = eye(opts_.ell) - [ zeros(1,opts_.ell-1) 0 ;
                                           eye( opts_.ell-1 ) zeros( opts_.ell-1 , 1 ) ];
            obj.Shm   = kron( prototype , eye(opts_.m) );
 
            obj.ShpInv = kron( prototype , eye(obj.pt) )\eye(obj.pt*opts_.ell); 
            obj.bcp    = [-eye(obj.pt); zeros((opts_.ell-1) * obj.pt, obj.pt)]; 
            

            
            
 
            %For the measurement constraints, we only care about the 
            %measurements that are part of our constraints so multiply
            %Gammahat and That by our C_l matrix which "picks" out
            %the measurements we want constrained
            if obj.numCon > 0 %if there are output constraints 
                Tc     = obj.C_l * obj.T;  
            else
                Tc     = [];
            end

            
            % We will add slack to yk constraints and NOT on uk mag constraints
            if ~obj.slackOn %no slack
                %Portion of A matrix for measurement constraints
                Abase  = Tc;
                %Portion of A matrix for U constraints  
                AextU   = eye(opts_.ell*opts_.m); 
                AextdU  = obj.Shm;   
            else %with slack
                %Portion of A matrix for measurement constraints
                Abase  = [Tc,  -eye(size(obj.S))  ]; 
                %Portion of A matrix for U constraints 
                AextU   = [ eye(opts_.ell*opts_.m) , zeros( opts_.ell*opts_.m , opts_.ell*obj.numCon ) ]; 
                %Portion of A matrix for deltaU constraints (takes Umpc and
                %makes deltaUmpc)
                AextdU   = [obj.Shm, zeros(obj.ell*obj.m, size(obj.S,2))];    
            end 
            %Portion of A matrix for introducing slack
            Aslack = [zeros( opts_.ell*obj.numCon , opts_.ell*opts_.m ),  -eye( opts_.ell*obj.numCon )];
             
            
            
            
            %Portion of b vector for measurement constraints
            alpha  = - kron(ones(obj.ell,1), obj.Dconstr); 
            %Portion of b vector for U constraints 
            bextU1 =  kron( ones(opts_.ell,1) , obj.u_max );
            bextU2 = -kron( ones(opts_.ell,1) , obj.u_min );
            %Portion of b vector for deltaU constraints
            if ~isempty(obj.Delu_max) 
                bext1 =  kron(ones(obj.ell,1), obj.Delu_max);
                bext2 = -kron(ones(obj.ell,1), obj.Delu_min); 
            end
            %Portion of b vector for introducing slack
            bslack = zeros( opts_.ell*obj.numCon , 1);

            
            if isempty(obj.Delu_max)
                AextdU = [];
                bext1  = [];
                bext2  = [];
            end
 
            %Initialize A b for Ax<=b constraints
            if obj.numCon == 0 %No output constraints
                obj.Aconstr = [ AextU; -AextU; AextdU; -AextdU ];
                obj.bconstr = [ bextU1; bextU2; bext1; bext2 ];  
            else
                if ~obj.slackOn %No slack
                    obj.Aconstr = [ Abase; AextU; -AextU; AextdU; -AextdU ];
                    obj.bconstr = [ alpha; bextU1; bextU2; bext1; bext2 ];  
                else %Slack active
                    obj.Aconstr = [ Abase; AextU; -AextU; AextdU; -AextdU; Aslack ];
                    obj.bconstr = [ alpha; bextU1; bextU2; bext1; bext2; bslack];  
                end
            end
             
            % Lagrange multipliers vector has dimensions of constraints.
            obj.lam     = 0*obj.bconstr; 
            
            if strcmp(opts_.optAlg,'admm') 
                % For admm slack has same dimensions as constraints so lam
                obj.slack = obj.lam; 
            end
              
            %Initialize MPC computed controls 
            obj.Umpc = zeros(opts_.m * opts_.ell,1);   
            obj.UC   = zeros(opts_.ell,1);
            %Set optimizer algorithm
            obj.optAlg = opts_.optAlg;
            obj.numIter = opts_.numIter;
                
            %flag for enabling explicit inverse
            obj.EGT_flag = opts_.EGT_flag; 
        end
        
        function updateUC(obj,UC)
            obj.UC = UC;
        end

        function  updateControl( obj, F_, G_, Y_ , U_ , r_ , specialFlag_ )
            %updateControl creates the new control input
             
            % Update Block matrices Fbar Gbar
            if ~obj.EGT_flag
                obj.updateFbarGbar( F_ , G_ , specialFlag_ );
            end
            
            % Update Command
            obj.updateCommand( r_ )  
  
            %quadratic MPC functions
            %state/measurement prediction
            obj.computeGammaT( Y_ , U_ , F_ , G_ , specialFlag_  );  

            % Update Tracking Error Integrator State
            r_temp = r_(:);
            obj.i = obj.i + ( obj.C_t*Y_(:,1) - r_temp(1:obj.pt) ) ;
            if specialFlag_
                [H,f] = obj.computeQPCostCoef( U_(1,1) );
                obj.computeQPConstraintCoef( U_(1,1) , specialFlag_ ); 
            else
                [H,f] = obj.computeQPCostCoef( U_(:,1) );
                obj.computeQPConstraintCoef( U_(:,1) ); 
            end

            %Quadratic Programming MPC solver
            switch obj.optAlg
                case 'quadprog'
                    [optSol, QPdiagnostics] = obj.computeQPresultQuadProg(H, f );
                case 'grad'
                    [optSol, QPdiagnostics] = obj.computeQPresultGrad(H, f );
                case 'admm'
                    [optSol, QPdiagnostics] = obj.computeQPresultADMM(H, f );
                otherwise
                    disp('Select quadprog, grad, or admm as solver.');
            end
            obj.extractQPControlandSlack( optSol , U_(:,1) , specialFlag_ );
 

        end
        
        %% Parameter Getters
        function [uOut,UOut] = getControl(obj)
            %getControl getter for the control input
            [uOut]  = obj.uOut;
            [UOut]  = obj.Umpc;
        end
    end
    
    
    methods (Access = private)                
        function updateCommand(obj, r_)
            %updateCommand Adds the newest command to the data buffers
          
            %update command 
            if length(r_) >  obj.pt 
                r_temp  = r_(:);
                obj.Rkl = r_temp(obj.pt+1:end);
            else
                obj.Rkl = kron(ones(obj.ell,1),r_);
            end
        end
        

        function updateFbarGbar(obj, F, G ,specialFlag)
            %updateFbarGbar adds newest F and G to Fbar and Gbar
 

            Frows = obj.p;
            Fcols = obj.p;
            Grows = obj.p;
            Gcols = obj.m;
 
            Fvec = zeros(Frows,Fcols*obj.nhat);
            Gvec = zeros(Grows,Gcols*(obj.nhat+1));
 
            for jj = 1:obj.nhat
                findices = (jj-1)*Fcols + 1 : jj*Fcols;
                gindices = (jj-1)*Gcols + 1 : jj*Gcols;
                Fvec(:,findices) = F(:,:,obj.nhat-jj+1 );
                Gvec(:,gindices) = G(:,:,obj.nhat-jj+2 ); 
            end
            gindices = obj.nhat*Gcols + 1 : (obj.nhat+1)*Gcols;
            Gvec(:,gindices) = G(:,:,1 );  

            for ii = 1:obj.ell
                obj.Fbar( (ii-1)*Frows+1 : ii*Frows , (ii-1)*Fcols+1 : (ii-1)*Fcols+1 + Fcols*(obj.nhat)-1 )    = Fvec;
                obj.Gbar( (ii-1)*Grows+1 : ii*Grows , (ii-1)*Gcols+1 : (ii-1)*Gcols+1 + Gcols*(obj.nhat+1)-1  ) = Gvec;
            end


 
            if specialFlag
 
                Grows = obj.p;
                Gcols = obj.m;
                Gvec = zeros(Grows,Gcols*(obj.nhat+1));
                for jj = 1:obj.nhat 
                    gindices = (jj-1)*Gcols + 1 : jj*Gcols; 
                    Gvec(:,gindices) = G(:,:,obj.nhat-jj+2 ); 
                end
                gindices = obj.nhat*Gcols + 1 : (obj.nhat+1)*Gcols;
                Gvec(:,gindices) = G(:,:,1 );  
                for ii = 1:obj.ell 
                    obj.Gbar( (ii-1)*Grows+1 : ii*Grows , (ii-1)*Gcols+1 : (ii-1)*Gcols+1 + Gcols*(obj.nhat+1)-1  ) = Gvec;
                end


                Gfine   = G(:,1,:);  
                Grows = obj.p;
                Gcols = 1;
                Gvec = zeros(Grows,Gcols*(obj.nhat+1));
                for jj = 1:obj.nhat 
                    gindices         = (jj-1)*Gcols + 1 : jj*Gcols; 
                    Gvec(:,gindices) = Gfine(:,:,obj.nhat-jj+2 ); 
                end
                gindices = obj.nhat*Gcols + 1 : (obj.nhat+1)*Gcols;
                Gvec(:,gindices) = Gfine(:,:,1 );  
                for ii = 1:obj.ell 
                    obj.Gbarf( (ii-1)*Grows+1 : ii*Grows , (ii-1)*Gcols+1 : (ii-1)*Gcols+1 + Gcols*(obj.nhat+1)-1  ) = Gvec;
                end


                Gcoarse = G(:,2,:);
                Grows = obj.p;
                Gcols = 1;
                Gvec = zeros(Grows,Gcols*(obj.nhat+1));
                for jj = 1:obj.nhat 
                    gindices         = (jj-1)*Gcols + 1 : jj*Gcols; 
                    Gvec(:,gindices) = Gcoarse(:,:,obj.nhat-jj+2 ); 
                end
                gindices = obj.nhat*Gcols + 1 : (obj.nhat+1)*Gcols;
                Gvec(:,gindices) = Gcoarse(:,:,1 );  
                for ii = 1:obj.ell 
                    obj.Gbarc( (ii-1)*Grows+1 : ii*Grows , (ii-1)*Gcols+1 : (ii-1)*Gcols+1 + Gcols*(obj.nhat+1)-1  ) = Gvec;
                end
            end
        end
           
        
        function computeGammaT(obj, Y, U, F, G , specialFlag )
            %computeGammaT computes Gamma and T matrices for response
            %prediction. 
            
            if specialFlag 
                Fcols = obj.p;  
    
                Fd = obj.Fbar(:,1:Fcols*obj.nhat);
                Fp = obj.Fbar(:,Fcols*obj.nhat+1:end);

                Gcols = obj.m; 
                Gd    = obj.Gbarc(:,1:Gcols*obj.nhat);  
                 
                %Update Gamma
                ytemp  = fliplr( Y(:,1:obj.nhat) );
                utemp  = fliplr( U(:,1:obj.nhat) );
                obj.Gamma = Fp\( Gd*utemp(:) - Fd*ytemp(:) );
                
                Gpc = obj.Gbarc(:,1*obj.nhat+1:end);
                obj.Tco = Fp\Gpc; 

                obj.Gamma = obj.Gamma + obj.Tco*obj.UC;


                %Update Block-Toeplitz Matrix 
                Gpf = obj.Gbarf(:,1*obj.nhat+1:end);
                obj.Tfi = Fp\Gpf;
                obj.T  = obj.Tfi;
            else 
                if obj.EGT_flag  
                    % Compute Gamma 
                    % Gamma0
                    Gamma0 = 0;
                    for jj = 1:obj.nhat
                        Gamma0 = Gamma0 + G(:,:,obj.nhat-jj+2)*U(:,obj.nhat-jj+1);
                    end
                    for jj = 1:obj.nhat-obj.mugap
                        Gamma0 = Gamma0 - F(:,:,obj.nhat-jj+1)*Y(:,obj.nhat-jj+1);
                    end
                    obj.Gamma( 1:obj.p , : ) = Gamma0;
                    % Gamma1 to Gammanhat-1
                    for ii = 1:obj.nhat-1
                        GammaI        = 0;
                        Gamma_indices = ii*obj.p + 1 : (ii+1)*obj.p;
                        for jj = 1:obj.nhat-ii
                            GammaI = GammaI + G(:,:,obj.nhat-jj+2)*U(:,obj.nhat-jj-ii+1);
                        end
                        for jj = 1:obj.nhat-max(obj.mugap,ii)
                            GammaI = GammaI - F(:,:,obj.nhat-jj+1)*Y(:,obj.nhat-jj-ii+1);
                        end        
                        for jj = 1:ii-obj.mugap
                            Gamma_indices2 = (jj-1)*obj.p + 1 : jj*obj.p;
                            GammaI = GammaI - F(:,:, ii - jj + 1 )*obj.Gamma( Gamma_indices2 , : ) ;
                        end
                        obj.Gamma( Gamma_indices , : ) = GammaI;
                    end
                    % Gammanhat to Gammaell-1
                    for ii = obj.nhat:obj.ell-1 
                        GammaI        = 0;
                        Gamma_indices = ii*obj.p + 1 : (ii+1)*obj.p;    
                        for jj = 1:obj.nhat-obj.mugap
                            Gamma_indices2 = (ii-obj.nhat+jj-1)*obj.p + 1 : (ii-obj.nhat+jj)*obj.p;
                            GammaI = GammaI - F(:,:,obj.nhat - jj + 1)*obj.Gamma( Gamma_indices2 , : ) ;
                        end
                        obj.Gamma( Gamma_indices , : ) = GammaI;
                    end
     
                    % Compute T
                    obj.T( 1:obj.p, 1:obj.m ) = G( : , : , 1 ); % H0
                    for ii = 1:obj.nhat 
                        TI = G( : , : , ii+1 ); % Gi
                        T_indices          = ii*obj.p + 1 : (ii+1)*obj.p;  
                        for jj = 1:ii-obj.mugap 
                            T_indices_2        = (jj-1)*obj.p + 1 : jj*obj.p;
                            TI = TI - F( : , : , ii-jj+1 )*obj.T( T_indices_2 , 1:obj.m );
                        end
                        obj.T( T_indices , 1:obj.m ) = TI;
                    end
                    for ii = obj.nhat+1:obj.ell-1 
                        TI = 0;
                        T_indices          = ii*obj.p + 1 : (ii+1)*obj.p;  
                        for jj = 1:obj.nhat-obj.mugap
                            T_indices_2        = (jj-1+ii-obj.nhat)*obj.p + 1 : (jj+ii-obj.nhat)*obj.p;
                            TI = TI - F( : , : , obj.nhat-jj+1 )*obj.T( T_indices_2 , 1:obj.m );
                        end
                        obj.T( T_indices , 1:obj.m ) = TI;
                    end
                    % Populate rest of the columns of T using the first column
                    for ii = 2:obj.ell
                        obj.T( (ii-1)*obj.p+1:end , (ii-1)*obj.m+1:ii*obj.m ) = obj.T( 1:end-(ii-1)*obj.p, 1:obj.m );
                    end 
                else
                    Fcols = obj.p;  
                    Gcols = obj.m;
        
                    Fd = obj.Fbar(:,1:Fcols*obj.nhat);
                    Fp = obj.Fbar(:,Fcols*obj.nhat+1:end);
                    Gd = obj.Gbar(:,1:Gcols*obj.nhat);
                    Gp = obj.Gbar(:,Gcols*obj.nhat+1:end);
                     
                    %Update Gamma
                    ytemp = fliplr( Y(:,1:obj.nhat) );
                    utemp = fliplr( U(:,1:obj.nhat) );
                    obj.Gamma = Fp\( Gd*utemp(:) - Fd*ytemp(:) );
                    
                    %Update Block-Toeplitz Matrix
                    obj.T = Fp\Gp; 
                end
            end
        end
        
        %% QP Methods
        function [H,f] = computeQPCostCoef(obj, u_ )
            %computeQPCostCoef Take the previous computed matrices
            %and combine them into a quadratic cost of the form
            %1/2*[Umpc;slack]'*H*Umpc + f'*[Umpc;slack] (1), where [Umpc;slack] is the
            %vector of variables to find in order to minimize the cost.
            %Umpc is the future control inputs and slack is a slack variable.
            %
            % 
            %deltaUmpc = Shm*Umpc-bc*u_k, Shm and bc are shift matrices 
            %
            %The slack variable 'slack' is (mpc horizon * num constraints x 1)
            %
            %u_ is the current control input
            %
            % The rows and columns corresponding to slack variables can be
            % removed by setting slackOn = 0.
             
            %For the cost function we only care about the tracked
            %measurements so multiply Gammahat and That by our C_tl matrix
            %which "picks" out the tracked measurements
            Tc     = obj.C_tl*obj.T;
            GammaC = obj.C_tl*obj.Gamma;

            %these variables come out of the conversion from cost (2) to
            %cost (1)
            alpha  = GammaC - obj.Rkl(:);
            z      = obj.bc*u_; 
            
            %compute the portion of H that is not multiplied by the slack
            %variable 'slack'   
            Hbar   = Tc'*obj.Q*Tc +  obj.Ru + obj.Shm'*obj.Rdu*obj.Shm;
            
            %compute the portion of f that is not multiplied by the slack  
            fbar   = ( Tc' * obj.Q * alpha +  obj.Shm' * obj.Rdu * z  ); 
 

            if det(obj.Qi) > 0
                % Tracking Error Integrator State
                TcI    = obj.ShpInv*obj.C_tl*obj.T;
                alphaI = obj.ShpInv*( GammaC - obj.Rkl(:) - obj.bcp*obj.i ); 

                %compute the portion of H that is not multiplied by the slack
                %variable 'slack'   
                Hbar   = Hbar + TcI'*obj.Qi*TcI ;
                
                %compute the portion of f that is not multiplied by the slack  
                fbar   = fbar + TcI' * obj.Qi * alphaI ; 
            end

            %create the full H variable
            %create the full f variable
            if obj.slackOn
                H      = blkdiag(Hbar, obj.S);
                f      = [fbar; zeros(size(obj.S,1), 1)];
            else
                H      = Hbar;
                f      = fbar;
            end
            
            H      = sparse( (H + H')/2 ); %ensure it is positive semidefinite
        end
        
        function computeQPConstraintCoef(obj, u_, varargin )
            %computeQPConstraintCoef Take the previous computed matrices
            %and update the constraints for the control and measurements in
            %quadratic program format:
            % A*[Xmpc;slack] <= b
            %
            %The constraints in the paper are in the format
            % C_l*Y1|k + D <= slack
            % Umpc_min <= Umpc <= Umpc_max
            % deltaUmpc_min <= deltaUmpc <= deltaUmpc_max
            % 0 <= slack
            %
            %
            %u_ is the current control input
            %
            % The rows and columns corresponding to slack variables can be
            % removed by setting slackOn = 0.
             
            if ~isempty(varargin)
                mtemp = 1;
            else
                mtemp = obj.m;
            end
            if ~isempty(obj.Delu_max)
                %Portion of b vector for deltaU constraints
                bext1 =  kron(ones(obj.ell,1), obj.Delu_max) - obj.bc * u_;
                bext2 = -kron(ones(obj.ell,1), obj.Delu_min) + obj.bc * u_;
            end

 

            %For the measurement constraints, we only care about the 
            %measurements that are part of our constraints so multiply
            %Gammahat and That by our C_l matrix which "picks" out
            %the measurements we want constrained
            if obj.numCon > 0 %if there are output constraints 
                Tc     = obj.C_l * obj.T; 
                GammaC = obj.C_l * obj.Gamma;
                %Portion of b vector for measurement constraints
                alpha  = -GammaC  - kron(ones(obj.ell,1), obj.Dconstr);
                
                % Only updates parts of Aconstr and bconstr that change with
                % data
                % Row arrangement for construction of A and b:
                %     [ for output constraints (if active) ;
                %       for contral mag max constraints;
                %       for control mag min constraints;
                %       for control move max constraints;
                %       for control move min constraints ];
                obj.Aconstr( 1:obj.ell*obj.numCon , 1:obj.ell*mtemp  ) = Tc;
                obj.bconstr( 1:obj.ell*obj.numCon , : ) = alpha; 
                if ~isempty(obj.Delu_max)
                    obj.bconstr( obj.ell*( obj.numCon + 2*mtemp ) + 1 : obj.ell*( obj.numCon + 4*mtemp ) , : ) = [ bext1  ; bext2 ]; 
                end
            else
                % Only updates parts of Aconstr and bconstr that change with
                % data
                % no output constraints  
                if ~isempty(obj.Delu_max)
                    obj.bconstr( obj.ell*(2*mtemp)+1 :  obj.ell*(4*mtemp) , : ) = [ bext1 ; bext2 ] ; 
                end
            end 
        end
        
        function [optSol, QPdiagnostics] = computeQPresultQuadProg(obj, H_, f_ )
            %computeQPresultQuadProg interface for the matlab QP solver.  
              
            %the slack component of the optimization initial guess to 0. 
            %The vector quadprog is solving for is [Xmpc;slack]
             
            options = optimoptions('quadprog','Display','off','MaxIterations',obj.numIter); 
                [optSol,Jopt,exitflag,output] = quadprog(H_, f_, ...
                                               obj.Aconstr,  obj.bconstr, [], [], ...
                                               [], [], obj.Umpc,...
                                               options);  

            if ( exitflag == -2 )
                disp('quadprog has returned infeasible problem token (exitflag = -2). Check Constraints.');
                optSol = [obj.Umpc; zeros(length(obj.S),1)]; %zero out the slack variable
            elseif ( exitflag == -3 )
                optSol = [obj.Umpc;zeros(length(obj.S),1)]; %zero out the slack variable
                disp('quadprog has returned unbounded problem token (exitflag = -3). Check Constraints.');
            else
            end
   
                QPdiagnostics.Jopt       = Jopt;
                QPdiagnostics.exitflag   = exitflag;
                QPdiagnostics.info       = output;
        end
%         
        function [optSol, QPdiagnostics] = computeQPresultGrad(obj, H_, f_ )
            %computeQPresultGradDescent compute optimizer using gradient
            %descent on dual problem




            
            exitflag = 0;
            if sprank(H_) == size(H_,1)
                Q1 = obj.Aconstr*(H_\obj.Aconstr');
                Q2 = obj.Aconstr*(H_\f_);
                df = Q1*obj.lam + Q2 + obj.bconstr; % gradient of Lagrangian
                alpha = 1/norm(Q1,'fro'); % step length



                %%% ADMM
                HADMM       = H_(1:end - obj.ell*obj.numCon,1:end- obj.ell*obj.numCon);
                Aconstraint = obj.Aconstr(1:end - obj.ell*obj.numCon,1:end- obj.ell*obj.numCon); 
                tempInverse = -(HADMM + alpha*transpose(Aconstraint)*Aconstraint)\eye(size(HADMM,1));
                JoptTrend   = zeros(2,obj.numIter); 
                lamADMM = obj.lam(1:end - obj.ell*obj.numCon);
                slackADMM = obj.slack(1:end - obj.ell*obj.numCon);
                if isempty(slackADMM)
                    slackADMM = lamADMM;
                end
                %%% ADMM 

                for k = 1:obj.numIter
                    obj.lam = max(obj.lam - alpha*df, 0); % steepest descent method
                    df = Q1*obj.lam + Q2 + obj.bconstr; % update Lagrangian
                    optSol = -H_\(f_ + obj.Aconstr'*obj.lam);
                    JoptTrend(2,k) = 0.5*optSol'*H_*optSol + f_'*optSol;





                    %%% ADMM
                    Aconstraint = obj.Aconstr(1:end - obj.ell*obj.numCon,1:end- obj.ell*obj.numCon);
                    bconstraint = obj.bconstr(1:end- obj.ell*obj.numCon,: );
                    HADMM       = H_(1:end - obj.ell*obj.numCon,1:end- obj.ell*obj.numCon);
                    fADMM       = f_(1:end- obj.ell*obj.numCon,: );
                    U         = tempInverse*( fADMM + alpha*transpose(Aconstraint)*( -slackADMM - bconstraint + lamADMM/alpha ) );
                    tempTerm  = Aconstraint*U - bconstraint;
                    slackADMM = tempTerm + lamADMM/alpha;
                    lamADMM   = lamADMM + alpha*( tempTerm - slackADMM  ); 

                    JoptTrend(1,k) = 0.5*U'*HADMM*U + fADMM'*U;
                    %%% ADMM

                end

                optSol = -H_\(f_ + obj.Aconstr'*obj.lam);
                Jopt = 0.5*optSol'*H_*optSol + f_'*optSol;
                if isnan(optSol)
                    optSol = zeros(size(H_,1),1);
                    Jopt = 0;
                    disp('grad has returned infeasible solution.');
                    exitflag = -1;
                end
            else
                warning('matrix H is singular');
                Jopt = 0;
                optSol = zeros(size(H_,1),1);
                exitflag = -2;
            end
            
        
            QPdiagnostics.Jopt       = Jopt;
            QPdiagnostics.exitflag   = exitflag;
            QPdiagnostics.JoptTrend  = JoptTrend;
        end

        function [optSol, QPdiagnostics] = computeQPresultADMM(obj, H_, f_ )
            %computeQPresultADNN compute optimizer using alternating
            %direction method of multipliers
            
            exitflag = 0;
            U        = obj.Umpc;
            if sprank(H_) == size(H_,1) 
                alpha = 1/norm(H_,'fro'); % step length  
                tempInverse = -(H_ + alpha*transpose(obj.Aconstr)*obj.Aconstr)\eye(size(H_,1));
                JoptTrend   = zeros(2,obj.numIter);



                %%% grad stuff
                Q1 = obj.Aconstr*(H_\obj.Aconstr');
                Q2 = obj.Aconstr*(H_\f_);
                df = Q1*obj.lam + Q2 + obj.bconstr; % gradient of Lagrangian
                alpha = 1/norm(Q1,'fro'); % step length 
                lamGD = obj.lam;
                %%% grad stuff

                for k = 1:obj.numIter
                    U         = tempInverse*( f_ + alpha*transpose(obj.Aconstr)*( -obj.slack - obj.bconstr + obj.lam/alpha ) );
                    tempTerm  = obj.Aconstr*U - obj.bconstr;
                    obj.slack = tempTerm + obj.lam/alpha;
                    obj.lam   = obj.lam + alpha*( tempTerm - obj.slack  ); 

                    JoptTrend(1,k) = 0.5*U'*H_*U + f_'*U;



                    lamGD = max(lamGD - alpha*df, 0); % steepest descent method
                    df = Q1*lamGD + Q2 + obj.bconstr; % update Lagrangian
                    optSolGD = -H_\(f_ + obj.Aconstr'*lamGD); 
                    JoptTrend(2,k) = 0.5*optSolGD'*H_*optSolGD + f_'*optSolGD;

                end

                optSol = U;
                Jopt   = 0.5*optSol'*H_*optSol + f_'*optSol;
                if isnan(optSol)
                    optSol = zeros(size(H_,1),1);
                    Jopt = 0;
                    disp('ADMM has returned infeasible solution.');
                    exitflag = -1;
                end
            else
                warning('matrix H is singular');
                Jopt = 0;
                optSol = zeros(size(H_,1),1);
                exitflag = -2;
            end
            
        
            QPdiagnostics.Jopt       = Jopt;
            QPdiagnostics.exitflag   = exitflag;
            QPdiagnostics.JoptTrend  = JoptTrend;
        end

%         
        function extractQPControlandSlack(obj, XoptSol_, u_, specialFlag_)
            %getControlandSlack get the control input and slack variable
            %from the solution of the QP. 

            if specialFlag_
                mtemp = 1; 
            else
                mtemp = obj.m;
            end

            U = XoptSol_(1:mtemp*obj.ell); 

            if obj.slackOn %Extract slack variable when slack is active
                obj.slack = XoptSol_(end-length(obj.S)+1:end);
            end
              
            % These constraints are useful only when solver is failing to
            % satisfy them 
            % Constrain entire optimizer not just first entry
            for ii = 1:obj.ell
                index = (ii-1)*mtemp + 1 : ii*mtemp; 
                uTemp = U(index,:);  
                if ii == 1
                    uLast  = u_;
                else
                    index2 = (ii-2)*mtemp + 1 : (ii-1)*mtemp; 
                    uLast  = U(index2,:);
                end

                for uu = 1:mtemp %loop over each control entry 
                    if uTemp(uu) > obj.u_max(uu) %magsat upper
                        uTemp(uu) = obj.u_max(uu);
                    elseif uTemp(uu) < obj.u_min(uu) %magsat lower
                        uTemp(uu) = obj.u_min(uu);
                    end

                    if ~isempty(obj.Delu_max)
                        move = uTemp - uLast; %compute control move
                        if move(uu) > obj.Delu_max(uu) %ratesat upper
                            uTemp(uu) = uLast(uu) + obj.Delu_max(uu);
                        elseif move(uu) < obj.Delu_min(uu) %ratesat lower
                            uTemp(uu) = uLast(uu) + obj.Delu_min(uu);
                        end 
                    end
                end
                U(index,:) = uTemp;  
            end 
            obj.uOut = U(1:mtemp);  

            obj.Umpc                            = XoptSol_; 
        end
        
    end
end


