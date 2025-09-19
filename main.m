clc; clear all;

% --- SETUP --- %

% Constants
consts.c = 3*10^8; % speed of light, [m/s]
consts.eps0 = (1/(36*pi))*10^(-9); % vacuum permittivity, [F/m]

% Adjustable parameters
L = 25; % length of cavity, [m]
D1 = 0.0254/2; % diameter of mirror 1, [m]
D2 = D1; % diameter of mirror 2, [m]
Rc1 = 50; % radius of curvature for mirror 1, [m]
Rc2 = Rc1; % radius of curvature for mirror 2, [m]
N = 2048; % number of mesh points along each dim of mesh grid
lambda = 1.064e-6; % laser wavelength, [m]
W = 2*D1; % domain half width, [m]
CFL = 0.0625;

% Grid
x = linspace(-W,W,N);
y = x;
dx = x(2) - x(1);
dy = dx;
%dz = CFL*4*k0*dx^2; % CFL-like condition, [m]
dz = L; % make each step a trip across the cavity, [m]
[X,Y] = meshgrid(x,y); % space domain

% Derived parameters
k0 = 2*pi/lambda; % freespace wavenumber, [m^-1]
cmask1 = (X.^2 + Y.^2 <= (D1/2)^2); % clipping mask mirror 1 (RHS)
cmask2 = (X.^2 + Y.^2 <= (D2/2)^2); % clipping mask mirror 2 (LHS)

% Set up mirror physical parameters for plotting
r1 = D1/2; % radius of mirror 1
r2 = D2/2; % radius of mirror 2
theta = linspace(0,2*pi,400);
x_circ1 = r1*cos(theta);
y_circ1 = r1*sin(theta);
x_circ2 = r2*cos(theta);
y_circ2 = r2*sin(theta);

% Input beam
w0 = 0.001; % input beam waist, [m]
E0 = exp(-(X.^2+Y.^2)/w0.^2); % input wave
E = E0;
I0 = 0.5*consts.eps0*consts.c*abs(E0).^2; % initial intensity, [W/m^2]

% Stability
g1 = 1 - L/Rc1; % stability parameter 1
g2 = 1 - L/Rc2; % stability paramter 2
g = g1*g2; % stability product, 0 < g < 1 for a stable cavity

% Set up frequency space
kx = (2*pi/(N*dx)) * (-N/2 : N/2-1);   % range from -pi/dx to +pi/dx
ky = kx;  % symmetric, since dx = dy
[KX, KY] = meshgrid(kx, ky);

% Free space transfer function of propagation 
H = exp(1i/(2*k0)*dz*(KX.^2+KY.^2));

% Mirror phase screens
rmask1 = exp(1i*k0*(X.^2+Y.^2)/(Rc1)); % reflection mask mirror 1 (RHS)
rmask2 = exp(1i*k0*(X.^2+Y.^2)/(Rc2)); % reflection mask mirror 2 (LHS)
%rmask1 = exp(1i*k0*(X.^2+Y.^2)/(2*Rc1));
%rmask2 = exp(1i*k0*(X.^2+Y.^2)/(2*Rc2));

% Simulation settings
save_interval = 1; % save frequency
step = 0; % initialize step 
Z_traveled = []; % initialize prop distance array
Z_position = []; % initialize intra-cavity position array
num_steps = round(L/dz); % number of steps needed for one trip across the cavity
Es = struct(); % initialize a struct for saving intermediate E fields

%%

% --- SIMULATION --- %

% Example of one round-trip

% --------------------------------(1)----------------------------------- %

% Beam enters cavity at +L/2 (RHS)

% Propagation from +L/2 (RHS) to -L/2 (LHS)
[step, Z_traveled, Z_position, E, Es] = R_L(step, Z_traveled, Z_position, E, Es, save_interval, num_steps, dz, L, H);

% Beam interacts with mirror 2 (LHS) mirror at -L/2
E = E.*cmask2.*rmask2;

plot_intensity(E, consts, Z_traveled, Z_position, step, X, Y, x_circ2, y_circ2); % take a look at the beam

% --------------------------------(2)----------------------------------- %

% Propagation from -L/2 back to +L/2
[step, Z_traveled, Z_position, E, Es] = L_R(step, Z_traveled, Z_position, E, Es, save_interval, num_steps, dz, L, H);

% Beam interacts with mirror 1 (RHS) mirror at +L/2
E = E.*cmask1.*rmask1;

plot_intensity(E, consts, Z_traveled, Z_position, step, X, Y, x_circ2, y_circ2); % take a look at the beam

%%

% Set how many round trips across the cavity you want to take.
% e.g. RHS --> LHS --> RHS = 1 round trip.

num_round_trips = 100;

% Prepare video writer
v = VideoWriter('cavity_propagation_3D_tight.mp4','MPEG-4');
v.FrameRate = 10; % adjust playback speed
open(v);

P0 = sum(sum(abs(E).^2)); % initial power

for i = 1:num_round_trips
    [step, Z_traveled, Z_position, E, Es] = R_L(step, Z_traveled, Z_position, E, Es, save_interval, num_steps, dz, L, H);
    E = E.*cmask2.*rmask2;
    [step, Z_traveled, Z_position, E, Es] = L_R(step, Z_traveled, Z_position, E, Es, save_interval, num_steps, dz, L, H);
    E = E.*cmask1.*rmask1;
    
    % Visualization
    I = 0.5*consts.eps0*consts.c*abs(E).^2;
    surf(X, Y, I/max(I(:)), 'LineStyle','none');
    ylabel('[m]')
    xlabel('[m]')
    title({
        sprintf('Laser Mode at Z = %.1f m', Z_traveled(step)), ...
        sprintf('Intra-Cavity Position = %.1f', Z_position(step))
    })
    hold on;
    plot3(x_circ2, y_circ2, zeros(size(x_circ2)), 'r-', 'LineWidth', 2); % plot mirror outline
    hold off;
    % Force white background
    set(gcf,'Color','w'); % figure background
    set(gca,'Color','w'); % axes background
    axis('square')
    axis tight;
    view(2)
    getframe();

    % Capture frame
    frame = getframe(gcf);
    writeVideo(v, frame);

    % Rescale 
    Pk(i) = max(max(abs(E).^2)); % peak intensity at this point
    E = E*sqrt(P0/sum(sum(abs(E).^2))); % rescale E
    BeamWidth(i) = sum(sum(sqrt(X.^2+Y.^2).*abs(E).^2))/P0;
end

close(v);
disp('Animation saved as cavity_propagation.mp4');

%% Post-processing

% Trackers used for both numeric and analytic
dz_post = 1; % little steps across the cavity  
zs = -L/2:dz_post:+L/2; % z positions, measured relative to the cavity center, z = 0, [m]
numSteps = length(zs); % number of steps to be taken for one trip across the cavity

% Numeric
H_post = exp(1i/(2*k0)*dz_post*(KX.^2+KY.^2)); % compute free space transfer function for a shorter distance
E_num = E; % Use the last saved E, numeric complex field at current plane
dz_post = 1; % little steps across the cavity    
Es_num = cell(1, numSteps); wz_num = zeros(1, numSteps); % initialize cell arrays to save numeric results
[Es_num, wz_num] = Save_Numeric_Results(E_num, Es_num, wz_num, numSteps, H_post, X, Y, dx, dy);  
fprintf('Numeric w0: %.3g mm \n', min(wz_num));

% Analytic 
wz_ana = zeros(1, numSteps); Iana = cell(1, numSteps); xbar = zeros(1, numSteps); ybar = zeros(1, numSteps);
sigma_x = zeros(1, numSteps); sigma_y = zeros(1, numSteps); sigma_r = zeros(1, numSteps);
zr = pi * w0_num^2 / lambda; % Rayleigh range from the numeric waist



%%
for k = 1:numSteps
    zrel = zs(k); % z relative to waist at cavity center
    wz_ana(k) = w0_num*sqrt(1 +(zrel/zr).^2); % analytical spot size

    % include curvature and optional Gouy phase for a correct complex field
    if zrel == 0
        Rz = Inf;
    else
        Rz = zrel * (1 + (zr^2 / zrel^2));  % R(z) = z * (1 + (zr^2 / z^2))
    end

    psi = atan(zrel / zr);   % Gouy phase (not used in intensity)
    % analytic complex field (amplitude + curvature) using 1/e^2 convention:
    % E(r,z) = (w0/w(z)) * exp(- r^2 / w(z)^2) * exp(-i * k0 * r^2 / (2 Rz) ) * exp(i * psi)
    r2 = X.^2 + Y.^2;
    if isfinite(Rz)
        curvature_phase = exp(-1i * k0 * r2 ./ (2 * Rz));
    else
        curvature_phase = 1;
    end

    Eana = (w0_num./wz_ana(k)).*exp(-r2./(wz_ana(k).^2)).* curvature_phase.*exp(1i * psi);
    Iana{k} = abs(Eana).^2; % analytic intensity

    % Moments
    Pana = sum(Iana{k}(:)) * dx * dy;
    xbar(k) = sum(sum(Iana{k} .* X)) * dx * dy / Pana;
    ybar(k) = sum(sum(Iana{k} .* Y)) * dx * dy / Pana;
    sigma_x(k) = sqrt(sum(sum(Iana{k} .* (X - xbar(k)).^2)) * dx * dy / Pana);
    sigma_y(k) = sqrt(sum(sum(Iana{k} .* (Y - ybar(k)).^2)) * dx * dy / Pana);
end

%% Graphing

figure;
plot(zs, +wz_ana, 'b', 'LineWidth', 1.5, 'DisplayName', 'Analytical'); hold on; 
plot(zs, -wz_ana, 'b', 'LineWidth', 1.5);
%plot(0, w0_num, 'ro', 'MarkerFaceColor','r', 'DisplayName','numeric w0');
plot(zs, +wz_num, 'r', 'LineWidth', 1.5, 'DisplayName', 'Numeric'); hold on; 
plot(zs, -wz_num, 'r', 'LineWidth', 1.5);
xlabel('z [m]');
ylabel('Transverse position [m]');
title('Beam cross section');
grid on;
axis tight;
