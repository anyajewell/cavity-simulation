clear all;
clc;
figure;

% Constants
c = 3e8; % [m/s]
eps0 = (1/(36*pi))*10^(-9); % vacuum permittivity, [F/m]

% Simulation parameters
Lx = 2; % Length of square transverse domain (one side), [m]
N = 511; % sampling number
dx = Lx/N; % step size 

% Initial conditions
Ld = 1064*1e-9; % Laser wavelength, [m]
w0 = 1e-1; % initial beam width, [m]
k0 = (2*pi)/Ld; % wavenumber
Zmax = 40e3; % destination of z, [m] 
Z0 = -Zmax; % starting location, [m]
pos = Z0;
t0 = 0; t(1) = t0; % starting time, [s]
v0 = 0; % starting transverse velocity, [m/s]
L = abs(Zmax - Z0); % cavity length, [m]
Nz = 1e2; % number of steps
dz = L/Nz; % step size, [m]
dt = dz/c; % time step, [s]
tmax = dt*Nz; % max time, [s]
Omega = 0.01; % rotational velocity, [rad/sec]
accel = 10000; % transverse acceleration, [m/s^2]
RTs = 1; % number of round trips to take

D1 = 0.5; % large size to reduce clipping, [m]
D2 = D1; % diameter of mirror 2, [m]
Rc1 = L*2; % radius of curvature for mirror 1, [m]
Rc2 = Rc1; % radius of curvature for mirror 2, [m]

% Set up mirror circles for plotting
r1 = D1/2; % radius of mirror 1
r2 = D2/2; % radius of mirror 2
theta = linspace(0,2*pi,400);
x_circ1 = r1*cos(theta); y_circ1 = r1*sin(theta); x_circ2 = r2*cos(theta); y_circ2 = r2*sin(theta);

% Calculations
zr = pi*w0^2/Ld; % Rayleigh range
wz = w0*sqrt(1+(L/zr).^2); % spot size, analytic solution

% Set up video
videoname = sprintf('Omega=%.3f_L=%.0fm.mp4', Omega, L);
[saveFolder, v] = Set_Up_Video(videoname);

% Set up simulation
[x, y, X, Y] = Create_Grid(N, Lx, dx);
Gau_ini = (1/(w0*pi*0.5))*exp(-(X.^2+Y.^2)./(w0^2));
z = linspace(Z0,Zmax,Nz); % z locations to evaluate field
t = linspace(t0,tmax,Nz); % timesteps
Gau = Gau_ini;
centerx(1) = trapz(trapz(X.*abs(Gau).^2))/trapz(trapz(abs(Gau).^2));
centery(1) = trapz(trapz(Y.*abs(Gau).^2))/trapz(trapz(abs(Gau).^2));
 
% Mirror masks: reflecting lens phase screens and clipping masks
rmask1 = exp(1i*k0*(X.^2+Y.^2)/(Rc1)); % reflection mask mirror 1 (RHS)
rmask2 = exp(1i*k0*(X.^2+Y.^2)/(Rc2)); % reflection mask mirror 2 (LHS)
cmask1 = (X.^2 + Y.^2 <= (D1/2)^2); % clipping mask mirror 1 (RHS)
cmask2 = (X.^2 + Y.^2 <= (D2/2)^2); % clipping mask mirror 2 (LHS)

%%
% Propagation

for a = 1:RTs
    % L --> R
    for i = 1:Nz
        pos = z(i); % current position
        if pos + dz > L/2 
            dz1 = Zmax-pos; 
            if dz1 == 0 % beam already at mirror
                break
            else % propagate by fractional step to L (?)
            [Gau, centerx, centery, v] = Prop(Gau, Omega, accel, dt, c, Ld, dx, dz1, x, y, z, X, Y, k0, v, centerx, centery, i); % propagate fractional step to mirror
            Gau = Gau.*cmask1.*rmask1; % clip and shape the beam
            dz2 = dz-dz1;
            [Gau, centerx, centery, v] = Prop(Gau, Omega, accel, dt, c, Ld, dx, dz2, x, y, z, X, Y, k0, v, centerx, centery, i); % propogate fractional step back
            end
        else % propagate by dz, as normal
            Omega = Omega; % opportunity to read/update Omega from other code
            accel = accel; % opportunity to read/update accel from other code
            [Gau, centerx, centery, v] = Prop(Gau, Omega, accel, dt, c, Ld, dx, dz, x, y, z, X, Y, k0, v, centerx, centery, i); % propagation loop
        end
    end

    % Interact with mirror (RHS)
    loss_frac = 1 - sum(abs(Gau.*cmask1).^2,'all') / sum(abs(Gau).^2,'all'); % loss this time
    Gau = Gau.*cmask1.*rmask1; % clip and shape the beam
    
    % Turn around
    Zmax = -Zmax;
    Z0 = -Zmax;
    pos = Z0;
    dz = -dz;
    z = linspace(Z0,Zmax,Nz);

    for i = 1:Nz
        pos = z(i); % current position
        if pos + dz < -L/2 
            dz1 = Zmax-pos; 
            if dz1 == 0 % beam already at mirror
                break
            else % propagate by fractional step to L (?)
            [Gau, centerx, centery, v] = Prop(Gau, Omega, accel, dt, c, Ld, dx, dz1, x, y, z, X, Y, k0, v, centerx, centery, i); % propagate fractional step to mirror
            Gau = Gau.*cmask2.*rmask2; % clip and shape the beam
            dz2 = dz-dz1;
            [Gau, centerx, centery, v] = Prop(Gau, Omega, accel, dt, c, Ld, dx, dz2, x, y, z, X, Y, k0, v, centerx, centery, i); % propogate fractional step back
            end
        else % propagate by dz, as normal
            Omega = Omega; % opportunity to read/update Omega from other code
            accel = accel; % opportunity to read/update accel from other code
            [Gau, centerx, centery, v] = Prop(Gau, Omega, accel, dt, c, Ld, dx, dz, x, y, z, X, Y, k0, v, centerx, centery, i); % propagation loop
        end
    end

    % Turn around
    Zmax = -Zmax;
    Z0 = -Zmax;
    pos = Z0;
    dz = -dz;
    z = linspace(Z0,Zmax,Nz);

end
    
% Rescale 
% Pk(i) = max(max(abs(E).^2)); % peak intensity at this point
% E = E*sqrt(P0/sum(sum(abs(E).^2))); % rescale E
% BeamWidth(i) = sum(sum(sqrt(X.^2+Y.^2).*abs(E).^2))/P0;

%%

% Analytic solution
x_rot = -Omega/c * (z.^2 + 0.5*L*z);
x_acc = (v0/c)*(z - Z0) + (accel/(2*c^2))*(z - Z0).^2;
x_ana = x_rot + x_acc; % analytic solution from ray optics

% Plotting
Plot_Results_vs_Analytic(centerx, x_ana, z) % plot results of propagation

close(v);





