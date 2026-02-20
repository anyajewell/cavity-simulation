clear all
clc

% 
UpperLowerLimits = [];
time = 50; 
Ts_f = 0.005; %Sample rate of fastest loop - also used as simulation rate
tau  = 40; %How many times PCAC fine is faster than PCAC coarse

Ts_c = Ts_f*tau;

stepsf = floor( time/Ts_f );
stepsc = floor( time/Ts_c );
loslength = 150e3; % distance between target and pursuer satellite in m
 

%
% Define the power function below. Basically replace the loss
%


%%% Loss model parameters (Given) %%%

D_ap = 0.4; % aperature diameter (m)
lambda = 1.064e-6; % wavelength (m)
L_range = loslength; % Distance between satellites (m)
alpha_L = 0.05; % dimensionless loss prameters (From Dr. Limbach)

theta0 = D_ap / L_range;
thetaD = lambda / D_ap;

% Calculate the minimum possible loss (according to the function)
L_min = alpha_L * (thetaD^2) / (theta0^2);

%% PCAC Coarse Settings
ellc  = 10;			% prediction/control horizon
nhat  = 4;			% estimator model order
mugap = 0;          % denominator gap
Ru    = 0;	        % U weight
Rdu   = 1000;	    % delta-U weight
Qi    = 0;	        % Integral Weight
P0    = 1e6;		% initial covariance   
[optsc] = PCAC_define_opts(nhat, mugap , ellc, P0, Ru , Rdu , Qi , -0.1 , 0.1 , -0.05, 0.05 , UpperLowerLimits ); 
% function [opts] = PCAC_define_opts(nhat,  mugap , ell, P0, Ru, Rdu , Qi , u_ll , u_ul , du_ll, du_ul , UpperLowerLimits ) % for my reference


%% PCAC Fine Settings
ellf  = 30;			% prediction/control horizon
nhat  = 4;			% estimator model order
mugap = 0;          % denominator gap
Ru    = 0;	        % U weight
Rdu   = 1;	        % delta-U weight
Qi    = 0;	        % Integral Weight
P0    = 1e6;		% initial covariance   
% [optsf] = PCAC_define_opts(nhat , mugap , ellf, P0, Ru , Rdu , Qi , -10 , 10 , -100, 100 , UpperLowerLimits ); 
[optsf] = PCAC_define_opts(nhat , mugap , ellf, P0, Ru , Rdu , Qi , -0.01 , 0.01 , -0.1*Ts_f, 0.1*Ts_f , UpperLowerLimits ); 

%% Preview Settings

% Add steps to account for use of preview
stepsf = stepsf + ellc*tau;
stepsc = stepsc + ellc;

% Error condition if inner loop horizon overshoots the total predicted
% outer loop control
ellIMax   =  tau*ellc; 
if ellf > ellIMax  
    disp(strcat('Choose inner horizon less than ', num2str(ellIMax,'%1.0f') )  );
end

 
%% Satellite with Gimbal State Space Model Dynamics %%

% Physical properties
Is = 10;        % Satellite moment of inertia [kg*m^2]
Ig = 0.5;       % Gimbal moment of inertia [kg*m^2]
bg = 1;      % Gimbal damping coefficient [N*m*s/rad]

% State matrix A (4x4)
% State vector: [theta_s; omega_s; theta_g_rel; omega_g_rel]
A = [0    1    0    0;                    % d(theta_s)/dt = omega_s
     0    0    0    bg/Is;                % d(omega_s)/dt coupling
     0    0    0    1;                    % d(theta_g)/dt = omega_g (relative)
     0    0    0    -(Is+Ig)*bg/(Is*Ig)]; % d(omega_g)/dt with damping

% Input matrix B (4x2)
% Input vector: [tau_s; tau_g]
B = [0           0;                       % theta_s not directly affected
     1/Is        -1/Is;                   % omega_s affected by both torques
     0           0;                       % theta_g not directly affected
     -1/Is       (Is+Ig)/(Is*Ig)];       % omega_g coupling

% Output matrix C - measuring satellite attitude and gimbal relative angle
C = [1  0  0  0;                         % y1 = satellite attitude (inertial)
     0  0  1  0];                        % y2 = gimbal angle (relative)

% Feedthrough matrix D (typically zero for this system)
D = [0  0;
     0  0];
 
sysSat = ss(A, B, C, D);
Gd  = c2d(sysSat,Ts_f); 
Ad = Gd.A;
Bd = Gd.B;
Cd = Gd.C;
Dd = Gd.D;


%% p trajectory (Target satellite motion) %%
k_array =  1:stepsf+1;  

% Linear trajectory
v_rel = 46; % m/s
omega_LOS = v_rel / loslength;
scale_factor = 1;

skc = scale_factor * omega_LOS * k_array*Ts_f;
% skc is the LOS bearing angle trajectory of the target relative to the pursuer (in radians)
% output function must be modified if skc is defined in radians vs meters
% (to ensure units are consistent)
% This setup is for LINEAR increase in radians, uncomment the plot below to
% see what it looks like

% % For debugging          
% figure
% plot(k_array*Ts_f, skc);
% xlabel('time'); ylabel('skc')

%% Memory allocation
x      = zeros(4,stepsf+1); 

qc     = zeros(1,stepsc+1); % Measured satellite attitude angle (where we are pointing) (rad)
pc     = zeros(1,stepsc+1); % Measured LOS angle of the target satellite (where we want to point) (rad)
df     = zeros(1,stepsf+1); % Saturated fine pointing error (rad)
df_us  = zeros(1,stepsf+1); % True (unsaturated) fine pointing error (rad)
 
uc     = zeros(1,stepsc+1); % Coarse loop input
uf     = zeros(1,stepsf+1); % Fine loop input

lambdac = ones(1,stepsc+1);
lambdaf = ones(1,stepsf+1);
thetac  = zeros(optsc.p,optsc.ltheta,stepsc+1);

L = zeros(1,stepsf+1);     % loss history

% need to articfically increase thetaf dimensions because IDing MIMO system
% while m = 1 is set in the opts file for other reasons
% Because the inner-loop pointing error depends on the outer-loop command, 
% the identified inner-loop model includes the effect of the coarse-loop dynamics via the previewed coarse command.
optsf.ltheta = (nhat-mugap)*1 + (nhat+1)*2; %length of theta vector 
thetaf  = zeros(optsf.p,optsf.ltheta,stepsf+1);

%Initialize PCAC
PCACc = PCAC.PcacClass('Linear_RLS', 'QP', optsc, 0,0);
PCACf = PCAC.PcacClass('Linear_RLS', 'QP', optsf, 0,1);


%% Simulation

jj = 1; %variable for coarse loop iteration  

% initialize data 
x(:,1)  = [0;0;0;0];
[pc(:,1),qc(:,1),df(:,1),df_us(:,1)] = output(x(:,1),skc(:,1),loslength);

% Initialize Loss 
delta_theta0 = df_us(:,1);
L(1) = (alpha_L * thetaD^2 + delta_theta0.^2) ./ (theta0^2 + delta_theta0.^2);

% Pre populate pc since preview is needed in the simulation
for ii = 1:stepsf
    if mod( ii , tau ) == 0 
        [pc(:,jj),~,~] = output(x(:,ii),skc(:,ii),loslength);
        jj             = jj + 1;
    end
end

jj = 1;
U_temp_inner = zeros( (ellf + 1)*tau , 1 );
logerror = 10;

% Loop iteration
for ii = 1:stepsf - tau*ellc
 
    % Plant propogation
    x(:,ii+1)  = Ad*x(:,ii) + Bd*[uc(:,jj); uf(:,ii)]; 
    % x = [theta_s; omega_s; theta_g; omega_g]
    % uc(:,jj) is the current coarse torque (which is the same for all
    % fine points between two coarse points)
    % uf(:,ii) is the current fine torque

    % Saturate movement of gimbal
    if abs( x(3,ii+1) ) > deg2rad(0.2) % arbitrary number, can change
        x(3,ii+1) = sign( x(3,ii+1) )*deg2rad(0.2) ;
        x(4,ii+1) = 0;
    end

    % The fine control loop needs a preview of coarse control commands
    % Extracting the preview of coarse commands for the fine loop
    % Start at ii-1 for fine pointing to have no offest (ii-1 when ii = 1
    % is 0)
    
    indices = rem( ii - 1 , tau ) + 2  : ellf +  rem( ii - 1 , tau ) + 1 ; 
    UC      = U_temp_inner( indices ); 
 
    %%%% COARSE LOOP %%%%%
    if mod( ii , tau ) == 0 
        [~,qc(:,jj),~]            = output(x(:,ii+1),skc(:,ii+1),loslength);
        [ uc(:,jj+1) , U_temp ]   = PCACc.oneStep( uc(:,jj), 1000*qc(:,jj) , 1000*pc(:,jj:jj+ellc)' , 0 , [] ); 
        thetac(:,:,jj+1)          = PCACc.idObj.getTheta;
        lambdac(:,jj+1)           = PCACc.idObj.getLambda; 
        U_temp_inner              = kron( [uc(:,jj);U_temp] , ones(tau,1) ); 
        logerror                  = log10(abs(rad2deg(pc(:,jj)-qc(:,jj))));
        jj                        = jj + 1;
    end

    %%%%% FINE LOOP %%%%%

    % Compute Loss to feed into PCAC for 'd'
    delta_theta = df_us(:,ii);
    Lk = sign(delta_theta) * ((alpha_L * thetaD^2 + delta_theta.^2) ./ (theta0^2 + delta_theta.^2)); % loss equation, signed gives better performance

    % The line below uses df(:,ii) 
    uf(:,ii+1)        = PCACf.oneStep( [uf(:,ii);uc(:,jj)] , Lk, 0, 1 , UC ); % Substitute df(:,ii) for Lk if you want to use loss
    % uf(:,ii+1)        = PCACf.oneStep( [uf(:,ii);uc(:,jj)] , POWER, REF_ES_ANG, 1 , UC ); % Substitute df(:,ii) for Lk if you want to use loss
    %uf(:,ii+1)        = PCACf.oneStep( [uf(:,ii);uc(:,jj)] , y, ref, 1 , UC ); 
    thetaf(:,:,ii+1)  = PCACf.idObj.getTheta;
    lambdaf(:,ii+1)   = PCACf.idObj.getLambda;   
    [~,~,df(:,ii+1),df_us(:,ii+1)] = output(x(:,ii+1),skc(:,ii+1),loslength);

    % Store loss
    delta_theta_next = df_us(:,ii+1);
    L(ii+1) = ((alpha_L * thetaD^2 + delta_theta_next.^2) ./ (theta0^2 + delta_theta_next.^2)); % Store NL Loss

    ii/( stepsf  - ellc*tau) 
end

%% Plots 

%%%%%% Fine Pointing Plots %%%%%%
figure(1)
clf
set(gcf, 'color', [1 1 1]) 

subplot(5,1,1)
% stairs( k_array*Ts_f , rad2deg(df_us) , 'm', 'linewidth',2)
hold on
stairs( k_array*Ts_f , rad2deg(df) , 'b', 'linewidth',2)
stairs( k_array*Ts_f , df*0 , 'r--' , 'linewidth', 2 )
hold off  
ylabel('$d$ (deg)','FontSize',12,'Interpreter','latex') 
grid on
xlim([0 time])
title('Fine Pointing','FontSize',16,'Interpreter','latex')
% ylim([-1 1])

subplot(5,1,2)
stairs( k_array*Ts_f , log10(abs(rad2deg(df))) , 'b', 'linewidth',2)  
grid on
ylabel('log$_{10}(|d|)$ (deg)','FontSize',12,'Interpreter','latex') 
xlim([0 time])


subplot(5,1,3)
stairs( k_array*Ts_f , uf(1,:) , 'b', 'linewidth',2)
hold on
stairs( [0 stepsf]*Ts_f , optsf.u_max(1)*[1 1] , 'k--' , 'linewidth', 1 )
stairs( [0 stepsf]*Ts_f , optsf.u_min(1)*[1 1] , 'k--' , 'linewidth', 1 )
hold off
ylabel('$\delta$','FontSize',12,'Interpreter','latex') 
grid on
xlim([0 time])


subplot(5,1,4)
stairs( k_array*Ts_f , lambdaf , 'b', 'linewidth',2) 
ylabel('$\lambda_k$','FontSize',12,'Interpreter','latex')
grid on
xlim([0 time])

subplot(5,1,5)
stairs( k_array*Ts_f , reshape(thetaf , optsf.ltheta , stepsf+1)' ,   'linewidth',2)
ylabel('$\theta_k$','FontSize',12,'Interpreter','latex')
xlabel('$t$ (sec)','FontSize',14,'Interpreter','latex')
grid on 
xlim([0 time])


%%%%%%% Coarse Pointing Plots %%%%%%%

j_array = 1:stepsc+1;

figure(2)
clf
set(gcf, 'color', [1 1 1]) 

% q plot
subplot(5,1,1)
% stairs( k_array*Ts_f , rad2deg(x(1,:)) , 'm', 'linewidth',2)
hold on
stairs( j_array*Ts_c , rad2deg(qc) , 'b', 'linewidth',2)
stairs( j_array*Ts_c , rad2deg(pc) , 'r--' , 'linewidth', 2 )
hold off  
grid on
ylabel('$q$ (deg)','FontSize',14,'Interpreter','latex') 
title('Coarse Pointing','FontSize',16,'Interpreter','latex')
xlim([0 time])

% Log Error
subplot(5,1,2)
stairs( j_array*Ts_c , log10(abs(rad2deg(qc-pc))) , 'b', 'linewidth',2)  
grid on
ylabel('log$_{10}(|p-q|)$ (deg)','FontSize',14,'Interpreter','latex') 
xlim([0 time])

% Coarse Input
subplot(5,1,3)
stairs( j_array*Ts_c , uc(1,:) , 'b', 'linewidth',2)
hold on
stairs( [0 stepsc]*Ts_c , optsc.u_max(1)*[1 1] , 'k--' , 'linewidth', 1 )
stairs( [0 stepsc]*Ts_c , optsc.u_min(1)*[1 1] , 'k--' , 'linewidth', 1 )
hold off
ylabel('$\tau$','FontSize',14,'Interpreter','latex') 
grid on
xlim([0 time])

% Coarse Forgetting
subplot(5,1,4)
stairs( j_array*Ts_c , lambdac , 'b', 'linewidth',2) 
ylabel('$\lambda_k$','FontSize',14,'Interpreter','latex')
grid on
xlim([0 time])

% Coarse Model
subplot(5,1,5)
stairs( j_array*Ts_c , reshape(thetac , optsc.ltheta , stepsc+1)' ,   'linewidth',2)
ylabel('$\theta_k$','FontSize',14,'Interpreter','latex')
xlabel('$t$ (sec)','FontSize',14,'Interpreter','latex')
grid on 
xlim([0 time])


% %%
% figure
% plot( j_array*Ts_c , abs(qc-pc) , 'b', 'linewidth',2)  ;
% grid on
% ylabel('zkc (deg)','FontSize',14,'Interpreter','latex');
% xlabel('time');
% title('qc-pc')
% xlim([0 time])

