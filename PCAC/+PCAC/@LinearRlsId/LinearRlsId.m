classdef LinearRlsId < handle
    %LinearId  

    properties (Access = private)
        ii %step
        m %number of inputs
        p %number of outputs
        nhat %identification order
        mugap %gap in identification model F coefficients
        theta %ID coefficients vector 
        lambda %forgetting factor
        ltheta %length of ID coefficients vector   
        Y %regressor for the measurements: num measurements x identification order
        U %regressor for the inputs: num inputs x identification order
        Z %regressor for the performance variable y-phi*theta (for VRF)

        % ell-step ID things
        ellId %horizon for ID
        Phi %regressor for the regressor nhat*(p+m)+m x ellId

        %ARX matrices
        Fbuffer %Buffer to store F matrices (p x p x nhat)
        Gbuffer %Buffer to store G matrices(p x m x nhat+1)
        
        %VRF parameters
        computeForgettingFactor; % VRF function handle
        VRF_parameters;          % Parameters for VRF function
        lambda_max;              % Upper bound for lambda

        %RLS parameters
        computeVRFupdate;        % RLS function handle
        RLS_parameters;          % Parameters for RLS function
        uncertaintyMatrix;       % P or R depending on RLS method used
        
    end
    
    properties (Access = public)
        getP;                    % Getter function for P
    end

    methods (Access = public)
        function obj = LinearRlsId(opts_,specialFlag_)
            %IDCLASS Construct an instance of this class and initialize its
            %properties
 
            
            %Initialize sizes
            obj.m      = opts_.m; %number of inputs
            if specialFlag_
                obj.m      = 2; %number of inputs
            end
            obj.p      = opts_.p; %number of outputs
            obj.nhat   = opts_.nhat; %Identification order
            obj.mugap  = opts_.mugap; %Identification F coefficient gap
            obj.lambda = 1; 

            
            %Initialize step
            obj.ii = 1;

            %Intialize ID coeffieicents and initial covaraiance to the 
            %correct size if only given a scalar
            if ( length(opts_.theta0) == 1 && length(opts_.P0) == 1)  
                obj.ltheta = (obj.nhat-obj.mugap)*obj.p + (obj.nhat+1)*obj.m;
                obj.theta  = opts_.theta0*ones(obj.p,obj.ltheta);
                P      = opts_.P0 * eye(obj.ltheta);
        
            elseif ( length(opts_.theta0) == 1 && length(opts_.P0) > 1) 
                obj.ltheta = (obj.nhat-obj.mugap)*obj.p + (obj.nhat+1)*obj.m;
                obj.theta  = opts_.theta0*ones(obj.p,obj.ltheta);
                P      = opts_.P0 * eye(obj.ltheta);
        
            elseif ( length(opts_.theta0) > 1 && length(opts_.P0) == 1)  
                obj.ltheta = (obj.nhat-obj.mugap)*obj.p + (obj.nhat+1)*obj.m;
                obj.theta  = opts_.theta0;
                P      = opts_.P0 * eye(obj.ltheta);
            else 
                obj.theta = opts_.theta0;
                P     = opts_.P0 * eye(obj.ltheta);
            end
        
            %Initialize VRF method and parameters
            if opts_.forgetting_method == "RmsWindows"
                obj.computeForgettingFactor = @computeForgettingFactor_RmsWindow;
                obj.VRF_parameters.eta = opts_.eta;
                obj.VRF_parameters.tau_n = opts_.tau_n;
                obj.VRF_parameters.tau_d = opts_.tau_d;

            elseif opts_.forgetting_method == "Ftest"
                obj.computeForgettingFactor = @computeForgettingFactor_Ftest;
                obj.VRF_parameters.eta = opts_.eta;
                obj.VRF_parameters.tau_n = opts_.tau_n;
                obj.VRF_parameters.tau_d = opts_.tau_d;
                obj.VRF_parameters.alpha = opts_.alpha;

            elseif opts_.forgetting_method == "RmsRegression"
                obj.computeForgettingFactor = @computeForgettingFactor_RmsRegression;
                obj.VRF_parameters.eta = opts_.eta;
                obj.VRF_parameters.regression_window = opts_.regression_window;
            else
                error("Invalid forgetting Method")
            end

            obj.lambda_max = opts_.lambda_max;

            % Initialize RLS method and parameters
            if opts_.RLS_method == "RLS"
                obj.getP = @obj.getPfromP;   % Covariance Getter
                obj.uncertaintyMatrix.P = P; % Covariance matrix for RLS
                obj.computeVRFupdate = @computeVRFupdate_RLS; % RLS Update Function

            elseif opts_.RLS_method == "ERRLS"
                obj.getP = @obj.getPfromR;   % Covariance Getter
                obj.uncertaintyMatrix.R = inv(P); % Information matrix for RLS
                obj.computeVRFupdate = @computeVRFupdate_ERRLS; % RLS Update Function
                obj.RLS_parameters.R_infinity = obj.prepareRinfinity(opts_.R_infinity,obj.ltheta); % R_infinity matrix
           
            elseif opts_.RLS_method == "CRRLS"
                obj.getP = @obj.getPfromP;   % Covariance Getter
                obj.uncertaintyMatrix.P = P; % Covariance matrix for RLS
                obj.computeVRFupdate = @computeVRFupdate_CRRLS; % RLS Update Function
                obj.RLS_parameters.R_infinity = obj.prepareRinfinity(opts_.R_infinity,obj.ltheta); % R_infinity matrix
                [V_,D_] = eig(obj.RLS_parameters.R_infinity); 
                D_ = diag(D_);
                obj.RLS_parameters.V = V_; % Eigenvectors of R_infinity
                obj.RLS_parameters.D = D_; % Eigenvalues of R_infinity
                
            else
                error("Invalid RLS Method")
            end

            %Initialize Y, U, Phi, and Z regressors
            obj.ellId = opts_.ellId; 
            maxbufferlength = max( obj.nhat+1 , obj.ellId );
            obj.Y   = zeros( obj.p , maxbufferlength );
            obj.U   = zeros( obj.m , obj.nhat+1 ); 
            obj.Phi = zeros( obj.ltheta , obj.ellId ); 


            if opts_.forgetting_method == "RmsWindows" || opts_.forgetting_method == "Ftest"
                obj.Z = zeros(1, opts_.tau_d+1);
            elseif opts_.forgetting_method == "RmsRegression"
                obj.Z = zeros(1, opts_.regression_window+1);
            else
                error("invalid Forgetting Method")
            end


            %Initialize -F and G buffers     
            obj.Fbuffer = zeros(obj.p, obj.p, opts_.nhat);
            obj.Gbuffer = zeros(obj.p, obj.m, opts_.nhat+1);
            
        end
        
        function  updateCoeff(obj, u_, y_)
            %updateCoeff updates the ID coefficients, this is what a user
            %would run 
            %update regressor
            phi           = obj.makeRegressor(u_, y_); 
            %compute the forgetting factor
            obj.lambda        = obj.computeForgettingFactor(obj,y_, phi); 
            %RLS update theta using regressor and forgetting factor
            obj.computeVRFupdate(obj); 
            %convert resulting theta to ARX format
            [obj.Fbuffer, obj.Gbuffer] = obj.makeFG(obj.theta, obj.nhat, obj.mugap, obj.m, obj.p);
            obj.ii = obj.ii + 1;
        end
        
        %% Parameter Getters
        function thetaOut = getTheta(obj)
            %getTheta outputs the current coefficient vector 
            thetaOut = obj.theta; 
        end
        
        function lambdaOut = getLambda(obj)
            %getTheta outputs the current coefficient vector  
            lambdaOut = obj.lambda;
        end
        
        function POut = getPfromP(obj)
            %getPfromP outputs the current RLS covariance 
            POut = obj.uncertaintyMatrix.P;
        end

        function POut = getPfromR(obj)
            %getPfromR outputs the current RLS covariance 
            POut = inv(obj.uncertaintyMatrix.R);
        end
        
        function [ Yout , Uout ] = getIObuffers(obj)
            %getIObuffers get current buffers of Y and U
            Yout = obj.Y(:,1:obj.nhat+1);
            Uout = obj.U;
        end
        
        function [ Fout , Gout ] = getFG(obj)  
            Fout = obj.Fbuffer;
            Gout = obj.Gbuffer;
        end
    end
    
    methods (Access = private)     
        %% Regressor Methods
        function phi = makeRegressor(obj, u_, y_)
            %makeRegressor creates the regressor phi for identification
            %consisting of past inputs and measurements and the newest input 
            %(u_) and output (y_).
            %updates phi regressor with phi
            
            %phi = [-y_{k-1}; ...; -y_{k-nhat}; u_{k}; ...; u_{k-nhat}]
            obj.updateBuff(u_, y_) %add the new control and measurement to the top of the regressor stack 
            yindices = (2+obj.mugap):obj.nhat+1;
            Ytmp = obj.Y(:,yindices);
            phi = [-Ytmp(:); obj.U(:)]; %create the full regressor matrix


            obj.Phi(:,2:end) = obj.Phi(:,1:end-1);
            obj.Phi(:,1)     = phi;
        end
        
        function updateBuff(obj, u_, y_)
            %updateBuff Adds the new control and measurements to the top of
            %the regressor stack
            obj.Y(:,2:end) = obj.Y(:,1:end-1);
            obj.Y(:,1)     = y_;
    
            obj.U(:,2:end) = obj.U(:,1:end-1);
            obj.U(:,1)     = u_; 
        end
        
        %% Forgetting Factor Methods 

        % RMS Window Methods
        function lambda = computeForgettingFactor_RmsWindow(obj, y_, phi_)
            %computeForgettingFactor compute the forgetting factor
            obj.updateZ(y_ - obj.theta * phi_ );
            gval   = obj.g_RmsWindow();
            us     = obj.unitStep(gval);
            lambda = 1/(1 + obj.VRF_parameters.eta * gval * us);

            % Saturate lambda to lambda max
            lambda = min(lambda,obj.lambda_max);
        end

        function gval = g_RmsWindow(obj)
            sd    = sum(obj.Z)/obj.VRF_parameters.tau_d; %denominator
            sn    = sum(obj.Z(1:obj.VRF_parameters.tau_n))/obj.VRF_parameters.tau_n; %numerator
            if ( sd == 0 )
                gval = 0;
            else
                gval  = sqrt(sn/sd)-1;
            end
        end

        % FTest Methods
        function lambda = computeForgettingFactor_Ftest(obj,y_,phi_)
            obj.updateZ(y_ - obj.theta * phi_);
            gval   = obj.g_Ftest();
            us     = obj.unitStep(gval);
            lambda = 1/(1 + obj.VRF_parameters.eta * gval * us);

            % Saturate lambda to lambda max
            lambda = min(lambda,obj.lambda_max);
        end

        function gval = g_Ftest(obj)%try: Multivariate analysis of variance
            %http://www.sci.rmutt.ac.th/stj/index.php/stj/article/view/46
            %https://en.wikipedia.org/wiki/Multivariate_analysis_of_variance
            %https://en.wikipedia.org/wiki/Hotelling%27s_T-squared_distribution
            %https://ncss-wpengine.netdna-ssl.com/wp-content/themes/ncss/pdf/Procedures/NCSS/Multivariate_Analysis_of_Variance-MANOVA.pdf
            %http://ibgwww.colorado.edu/~carey/p7291dir/handouts/manova1.pdf
            %https://support.sas.com/documentation/cdl/en/statug/63033/HTML/default/viewer.htm#statug_introreg_sect012.htm
            %sd    = sum(obj.Z)/obj.tau_d; %denominator
            %sn    = sum(obj.Z(1:obj.tau_n+1))/obj.tau_n; %numerator
            sdtemp = obj.Z-mean(obj.Z,2);
            %sd = sum(vecnorm(sdtemp,2,1))/obj.tau_d;
            sd = sdtemp*sdtemp'/obj.VRF_parameters.tau_d;
            sntemp = obj.Z(:,1:obj.VRF_parameters.tau_n+1)-mean(obj.Z(:,1:obj.VRF_parameters.tau_n+1),2);
            %sn = sum(vecnorm(sntemp,2,1))/obj.tau_n;
            sn = sntemp*sntemp'/obj.VRF_parameters.tau_n;
            
            % PARAMETER:
            alpha = obj.VRF_parameters.alpha;

            if (det(sd)^(1/obj.p) < 1e-10) || obj.Z(end) == 0; 100*eps; %( sd == 0 )
                gval = 0;
            else
                %Lawley - Hotelling Trace
                h = obj.VRF_parameters.tau_n; 
                e = obj.VRF_parameters.tau_d; %dof
                a = obj.p*h;
                B = ((e+h-obj.p-1)*(e-1))/((e-obj.p-3)*(e-obj.p));
                b = 4+(a+2)/(B-1);
                c = a*(b-2)/(b*(e-obj.p-1));
                H = (obj.VRF_parameters.tau_n)*sn; %sum of squares formulation of variance
                E = (obj.VRF_parameters.tau_d)*sd; %sum of squares
                T = trace(H*E^-1); 
                F = T/c;
                
                testvar = F;
                
                %pool = ((obj.tau_n)*sn+(obj.tau_d)*sd)/((obj.tau_n)+(obj.tau_d));
                %testvar = ((obj.tau_n+1)*(obj.tau_d+1)/(obj.tau_d+1+obj.tau_n+1))*(mean(obj.Z(:,1:obj.tau_n+1),2)-mean(obj.Z,2))'*pool^-1*(mean(obj.Z(:,1:obj.tau_n+1),2)-mean(obj.Z,2));
                %testvar = ((obj.tau_n+1+obj.tau_d+1-1-obj.p)/((obj.tau_n+1+obj.tau_d+1-2)*obj.p))*testvar;
                
                if testvar < 1 %F-Test, %only care about if the numerator is larger?
                    gval  = -10; %sd/sn - finv(0.999999,obj.tau_d,obj.tau_n);
                else
                    gval  = sqrt(testvar) - sqrt(finv(1-alpha,a,b)); %variance was computed with obj.tau_n/d + 1 elements
                end
                if gval > 0
                    %gval = sqrt(gval);
                end
            end
        end



        % RMS Regression Methods
        function lambda = computeForgettingFactor_RmsRegression(obj,y_,phi_)
                
            % Update error 
            obj.updateZ(y_ - obj.theta * phi_);

            % Z(1) is most recent value, Z(end) is first value
            B = sqrt(obj.Z(1:obj.VRF_parameters.regression_window));  
            
            % B(1) is oldest data, abs_error(end) is newest
            B = fliplr(B);
            
            % Linear regression over abs_error
            A = [ones(size(B))' [1:obj.VRF_parameters.regression_window]'];
            X = A\B';

            slope = X(2);
            us     = obj.unitStep(slope);
            
            lambda = 1/(1 + obj.VRF_parameters.eta * slope * us);

            % Saturate lambda to lambda max
            lambda = min(lambda,obj.lambda_max);
        end

        % Common Methods
        function updateZ(obj, z_)
            obj.Z     = [norm(z_)^2 obj.Z(1:(end-1))];
        end

        function us = unitStep(~, x_) 
            if x_ < 0
                us = 0;
            else
                us = 1; 
            end
        end


        %% VRF Update methods
        function computeVRFupdate_RLS(obj)
            phi_ = obj.Phi';
            y_   = obj.Y(:,1:obj.ellId)'; 

            L     = obj.uncertaintyMatrix.P/obj.lambda;
            M     = phi_ * L * phi_';
            N     = L*phi_';
            obj.uncertaintyMatrix.P = L - N / (eye(size(M)) + M) * N';    
            obj.theta = obj.theta + ( obj.uncertaintyMatrix.P * phi_' * ( y_ - phi_ * ( obj.theta') ) )';    
        end

        function computeVRFupdate_ERRLS(obj) % Exponential Resetting RLS
            phi_ = obj.Phi';
            y_ = obj.Y(:,1:obj.ellId)';

            obj.uncertaintyMatrix.R = obj.lambda*obj.uncertaintyMatrix.R + (1-obj.lambda)*obj.RLS_parameters.R_infinity + phi_'*phi_;
            obj.theta = obj.theta + ( obj.uncertaintyMatrix.R \phi_' * (y_ - phi_ * obj.theta') )';   
        end

        function computeVRFupdate_CRRLS(obj) % Cyclic Resetting RLS
            phi_ = obj.Phi';
            y_ = obj.Y(:,1:obj.ellId)';

            k_cyclic = mod(obj.ii-1,obj.ltheta);

            vbar = sqrt(  (1-obj.lambda^obj.ltheta)*obj.RLS_parameters.D(k_cyclic+1)/(obj.lambda^(obj.ltheta-k_cyclic-1))  )*obj.RLS_parameters.V(:,k_cyclic+1);
            phibar = [phi_; vbar'];

            obj.uncertaintyMatrix.P = (  obj.uncertaintyMatrix.P - (obj.uncertaintyMatrix.P*phibar')/(obj.lambda*eye(obj.ellId+1) + phibar*obj.uncertaintyMatrix.P*phibar')*(phibar*obj.uncertaintyMatrix.P)  )/obj.lambda;
            obj.theta = obj.theta + ( obj.uncertaintyMatrix.P * phi_' * ( y_ - phi_ * ( obj.theta') ) )';   
        end
        
    end
    
    methods (Static)
        %Useful functions you can call without making an object
        [F, G] = makeFG(theta, nhat, mugap, m, p)

        function R_infinity = prepareRinfinity(R_infinity_,ltheta)
            if isscalar(R_infinity_)
                R_infinity = R_infinity_*eye(ltheta);
            elseif isequal(size(R_infinity_), [ltheta,ltheta])
                R_infinity = R_infinity_;
            else
                error("Invalid size for opts_.R_infinity");
            end
        end
    end
end

