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
Omega = 0; % angular rotation between spacecraft frame and inertial, geocentric frame, [rad/s]

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
w0 = 0.002; % input beam waist, [m]
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

% Propagation operators
H = exp(1i/(2*k0)*dz*(KX.^2+KY.^2)); % free space transfer function of propagation 
R = @(z1, z2) exp(-i*KX*Omega/consts.c*1/2*(z2+z1)*(z2-z1)*dz); % rotation operator

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
v = VideoWriter('cavity_propagation_rotating.mp4','MPEG-4');
v.FrameRate = 10; % adjust playback speed
open(v);

P0 = sum(sum(abs(E).^2)); % initial power

for i = 1:num_round_trips
    [step, Z_traveled, Z_position, E, Es] = R_L(step, Z_traveled, Z_position, E, Es, save_interval, num_steps, dz, L, H, R);
    E = E.*cmask2.*rmask2;
    [step, Z_traveled, Z_position, E, Es] = L_R(step, Z_traveled, Z_position, E, Es, save_interval, num_steps, dz, L, H, R);
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
    E = E*sqrt(P0/sum(sum(abs(E).^2))); % rescale E to have no diffraction losses
    BeamWidth(i) = sum(sum(sqrt(X.^2+Y.^2).*abs(E).^2))/P0;
end

close(v);

%% Post-processing

% Trackers used for both numeric and analytic
dz_post = .5; % little steps across the cavity  
zs = -L/2:dz_post:+L/2; % z positions, measured relative to the cavity center, z = 0, [m]
numSteps = length(zs); % number of steps to be taken for one trip across the cavity

% Numeric
H_post = exp(1i/(2*k0)*dz_post*(KX.^2+KY.^2)); % compute free space transfer function for a shorter distance
E0_num = E; % Use the last saved E, numeric complex field at current plane
dz_post = 1; % little steps across the cavity    
Es_num = cell(1, numSteps); wz_num = zeros(1, numSteps); % initialize cell arrays to save numeric results
[Es_num, wz_num] = Save_Numeric_Results(E0_num, Es_num, wz_num, numSteps, H_post, X, Y, dx, dy);  
I_num = cellfun(@(E) abs(E).^2, Es_num, 'UniformOutput', false); % numeric intensity
I0_num = abs(E0_num).^2; 
w0_num = min(wz_num); % numeric waist
fprintf('Numeric w0: %.3g mm \n', w0_num*10^3);
I_num = I0_num*exp(-2*(X.^2+Y.^2)/w0_num^2); % intensity at numeric 1/e^2 waist

% Analytic 
R = Rc1; % set which mirror the beam is arriving at
w0_ana = sqrt(lambda*L/pi * sqrt(g1*g2*(1-g1*g2)/(g1+g2-2*g1*g2)^2)); % analytic beam waist from cavity geometry
fprintf('Analytic w0: %.3g mm \n', w0_ana*10^3);
zr = pi * w0_ana^2 / lambda; % Rayleigh range from the numeric waist
wz_ana = w0_ana * sqrt(1+(zs./zr).^2); % analytic spot size along zs
q0 = 1i*(pi*w0_ana^2)/lambda; % initial complex radius, evaluated at the center of the cavity, i.e. at the beam waist
qz = R - (pi*wz_ana.^2)/(1i*lambda); % complex radius of curvature

% for i = 1:length(qz)
%     qzi = qz(i);
%     %E_ana{i} = 1/qzi * exp(-1i*k0 * (X.^2+Y.^2)/(2*qzi)); % analytical field at point z
%     %I_ana{i} = cellfun(@(E) abs(E).^2, E_ana, 'UniformOutput', false); % analytic intensity
% end

%% Graphing

figure;
plot(zs, +wz_ana, 'b', 'LineWidth', 1.5, 'DisplayName', 'Analytic'); hold on; 
plot(zs, -wz_ana, 'b', 'LineWidth', 1.5, 'HandleVisibility', 'off');
plot(zs, +wz_num, 'r', 'LineWidth', 1.5, 'DisplayName', 'Numeric'); hold on; 
plot(zs, -wz_num, 'r', 'LineWidth', 1.5, 'HandleVisibility', 'off');
xlabel('z [m]');
ylabel('Transverse position [m]');
title('Beam cross section');
grid on;
legend()
axis tight;
