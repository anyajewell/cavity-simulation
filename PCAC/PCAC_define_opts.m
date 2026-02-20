function [opts] = PCAC_define_opts(nhat,  mugap , ell, P0, Ru, Rdu , Qi , u_ll , u_ul , du_ll, du_ul ,  UpperLowerLimits )
% PCAC_define_opts: Define PCAC parameters for illustrative MIMO example
%
% Inputs:
% nhat             -- Model order
% mugap            -- denominator gap in SISO, number of F matrices zeroed for MIMO
% ell              -- Horizon
% P0               -- Recursive Least Squares initial covariance P0*eye(ltheta)
% Ru               -- Control weight
% Rdu              -- Control move weight
% Qi               -- Integration weight
% u_ll             -- Control lower limit
% u_ul             -- Control upper limit
% du_ll            -- Control-move lower limit
% du_ul            -- Control-move upper limit
% UpperLowerLimits -- Output constraints
% 
% Outputs:
% opts             -- struct containing PCAC parameters
%
% Usage: [opts] = PCAC_define_opts()
%
% Author: Syed Aseem Ul Islam and Jacob Vander Schaaf 
% Last Modified: March 1, 2024
%
% Required files: None. 
if mugap > nhat
    error('mugap must be smaller than or equal to nhat')
end
if mugap == nhat
    warning('PCAC will implement an FIR model')
end
 
 

%% Structural Parameters
p = 1; % number of outputs 
m = 1; % number of inputs 
 
 

%% Estimator Parameters
theta0 = 0.01;           % initial parameter estimate 
ltheta = (nhat-mugap)*p + (nhat+1)*m; %length of theta vector 

ellId = 1; %Horizon used for RLS ID

 
EGT_flag = 0; %Flag to enable explicit computation of Gamma and T matrices
 

% RLS
% RLS: Recursive Least Squares
% ERRLS: Exponential Resetting RLS
% CRRLS: Cyclic Resetting RLS 
RLS_method = "RLS";           % ERRLS or CRRLS or RLS
% RLS_method = "CRRLS";        % ERRLS or CRRLS or RLS
% RLS_method = "ERRLS";        % ERRLS or CRRLS or RLS

% For ERRLS or CRRLS
R_infinity = 0.000001;


lambda_max = 1; % Lambda upper bound

% VRF
forgetting_method = "RmsWindows";    % RmsWindows, Ftest, or RmsRegression
% forgetting_method = "Ftest";         % RmsWindows, Ftest, or RmsRegression
% forgetting_method = "RmsRegression"; % RmsWindows, Ftest, or RmsRegression

% Forgetting parameter
eta = 0.01;

% For RmsWindows or Ftest
tau_n = 40;
tau_d = 200;

% For Ftest
alpha = 0.01;

% For RmsRegression
regression_window = 5;



%% Output Parameters
C_t = eye(p);   % select tracking outputs    % pt x p % set to the identity if tracking all feedback signals
 

% C_c = 1*eye(p); % select constrained outputs % pc x p % currently, just set pc = p and zero out. 
% C_c = [1,2;0,1];
C_c = []; %testing no output constraints



if isempty(C_c) %no output constraints
    C = [];
    D = [];
else
    % output constraints C * yc + D <= 0 
    C = 0*eye(p); 
    C = kron( eye(size(C_c,1)) ,[1;-1] ); %lets us set upper and lower bounds simultaneously in the C*yc + D <= 0 form

    D = 0*ones(p,1); % C in nc x pc, D in nc x 1 
    D = kron(  eye(size(C_c,1)),[-1 0;0 1] )*UpperLowerLimits;
end

pt = size(C_t,1);
nc = length(D);

%% Optimization Paramters  
Qhat = 1 * eye( (ell-1) * pt );      % QP cost-to-go weight     % (ell-1)*pt x (ell-1)*pt
Phat = 1 * eye( pt );                % QP terminal value weight %  pt x pt
Rdu  = Rdu * eye( ell * m );         % QP delta-U weight        %  ell*m x ell*m   
Ru   = Ru * eye( ell * m );          % QP U weight              %  ell*m x ell*m   
S    = 1 * eye( ell * nc );          % QP slack weight          %  ell*nc x ell*nc
% S    = 0 * eye( ell * nc ); % QP slack weight set to 0 % ell*nc x ell*nc


u_min    = u_ll*ones(m,1); % m x 1
u_max    = u_ul*ones(m,1); % m x 1
if ~( isempty(du_ll) || isempty(du_ul) )
    Delu_min = du_ll*ones(m,1); % input constraints % m x 1 % pm 1
    Delu_max = du_ul*ones(m,1); % m x 1
else
    Delu_min = []; % input constraints % m x 1 % pm 1
    Delu_max = []; % m x 1
end

%optimization algorithm. Use Matlab's quadprog or a projected gradient
%algorithm
optAlg = 'quadprog'; %quadprog or grad or admm
numIter = 200; %number of iterations. max iterations for quadprog (default for quadprog is 200)


%% create opts struct

% Universal options
opts = struct('p',        p,...
              'm',        m,...
              'nhat',     nhat,...
              'mugap',    mugap,...
              'ell',      ell,...
              'theta0',   theta0,...
              'ltheta',   ltheta,...
              'P0',       P0,... 
              'ellId',    ellId,...
              'eta',      eta,...
              'C_t',      C_t,...
              'C_c',      C_c,...
              'C',        C,...
              'Dconstr',  D,...
              'Qhat',     Qhat,...
              'Phat',     Phat,...
              'Ru',       Ru,...
              'Rdu',      Rdu,...
              'Qi',       Qi,...
              'S',        S,...
              'Delu_min', Delu_min,... 
              'Delu_max', Delu_max,... 
              'u_min',    u_min,...    
              'u_max',    u_max,...
              'optAlg',   optAlg,...
              'numIter',  numIter,...
              'forgetting_method', forgetting_method,...
              'RLS_method', RLS_method,...
              'lambda_max', lambda_max,...
              'EGT_flag',EGT_flag);
           
if RLS_method == "ERRLS" || RLS_method == "CRRLS"
    opts.R_infinity = R_infinity;
end

if forgetting_method == "RmsWindows" || forgetting_method == "Ftest"
    opts.tau_n = tau_n;
    opts.tau_d = tau_d;
    if forgetting_method == "Ftest"
        opts.alpha = alpha;
    end
elseif forgetting_method == "RmsRegression"
    opts.regression_window = regression_window;
end

end