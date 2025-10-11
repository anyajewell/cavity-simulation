clc; clear all;

% --- SETUP --- %

% Constants
consts.c = 3*10^8; % speed of light, [m/s]
consts.eps0 = (1/(36*pi))*10^(-9); % vacuum permittivity, [F/m]

% Adjustable parameters
L = 500; % length of cavity, [m]
D1 = 0.0254/2; % diameter of mirror 1, [m]
D2 = D1; % diameter of mirror 2, [m]
Rc1 = 2000; % radius of curvature for mirror 1, [m]
Rc2 = Rc1; % radius of curvature for mirror 2, [m]
N = 2048; % number of mesh points along each dim of mesh grid
lambda = 1.064e-6; % laser wavelength, [m]
W = 4*D1; % domain half width, [m]
CFL = 0.0625; % CFL number
Omega = 0; % relative rotation of spacecraft frame to inertial geocentric frame, [rad/s]

% Grid
k0 = 2*pi/lambda; % freespace wavenumber, [m^-1]
x = linspace(-W,W,N);
y = x;
dx = x(2) - x(1);
dy = dx;
dz = CFL*4*k0*dx^2; % CFL-like condition, [m]
dz = L; % make each step a trip across the cavity, [m]
%dz = 1;
[X,Y] = meshgrid(x,y); % space domain

% Set up mirror physical parameters for plotting
r1 = D1/2; % radius of mirror 1
r2 = D2/2; % radius of mirror 2
theta = linspace(0,2*pi,400);
x_circ1 = r1*cos(theta); y_circ1 = r1*sin(theta); x_circ2 = r2*cos(theta); y_circ2 = r2*sin(theta);

% Input beam
w0 = 0.001; % input beam waist, [m]
zr = pi*w0^2/lambda;
E0 = exp(-(X.^2+Y.^2)/w0.^2); % input wave
E = E0;
I0 = 0.5*consts.eps0*consts.c*abs(E0).^2; % initial intensity, [W/m^2]

% Stability
g1 = 1 - L/Rc1; % stability parameter 1
g2 = 1 - L/Rc2; % stability paramter 2
g = g1*g2; % stability product, 0 < g < 1 for a stable cavity

% Set up frequency space
kx = (2*pi/(N*dx)) * (-N/2 : N/2-1); % range from -pi/dx to +pi/dx
ky = kx;  % symmetric, since dx = dy
[KX, KY] = meshgrid(kx, ky);

% Propagation operators, to be used in frequency space
H = exp(1i/(2*k0)*dz*(KX.^2+KY.^2)); % free space transfer function of propagation
R = @(z1, z2) exp(-i*KX*Omega/consts.c*1/2*(z2+z1)*dz); % rotation operator
% 2) symmetric split-step propagation (single dz step example)
% Precompute H on centered K-grid (matching fftshift convention)
kx = (2*pi/(N*dx)) * (-N/2 : N/2-1);
[KX, KY] = meshgrid(kx, kx);
H = exp(1i/(2*k0)*dz*(KX.^2 + KY.^2));

% Propagation masks: mirror phase screens, clipping masks, and tilting
% mask, to be used in real space
rmask1 = exp(1i*k0*(X.^2+Y.^2)/(Rc1)); % reflection mask mirror 1 (RHS)
rmask2 = exp(1i*k0*(X.^2+Y.^2)/(Rc2)); % reflection mask mirror 2 (LHS)
cmask1 = (X.^2 + Y.^2 <= (D1/2)^2); % clipping mask mirror 1 (RHS)
cmask2 = (X.^2 + Y.^2 <= (D2/2)^2); % clipping mask mirror 2 (LHS)
%cmask1 = 1; cmask2 = cmask1; % turn off clipping
tmask = exp(k0*Omega*X./(1i*consts.c)*dz); % tilting mask, derived by me
%tmask = 1; % turn off t mask

% Set up an absorbing mask
ra = sqrt(X.^2 + Y.^2); ra_norm = ra / W; amask = ones(size(X));
edge_start = 0.8; % start absorbing after 80% of domain
mask_zone = ra_norm > edge_start;
amask(mask_zone) = cos((pi/2) * (ra_norm(mask_zone) - edge_start) / (1 - edge_start)).^2;
amask(ra_norm >= 1) = 0; % fully absorb at edge

% Simulation settings
save_interval = 1; % save frequency
step = 0; % initialize step 
Z_traveled = []; % initialize prop distance array
Z_position = []; % initialize intra-cavity position array
num_steps = round(L/dz); % number of steps needed for one trip across the cavity
Es = struct(); % initialize a struct for saving intermediate E fields

% Prepare video writer
todayStr = datestr(now, 'yyyy-mm-dd');
saveFolder = fullfile('C:\Users\Anya Jewell\Documents\MATLAB\ORACLE\Results', todayStr);

if ~exist(saveFolder, 'dir')
    mkdir(saveFolder);
end

fileName = sprintf('Omega=%.3f_L=%.0fkm.mp4', Omega, L);
filePath = fullfile(saveFolder, fileName);
v = VideoWriter(filePath, 'MPEG-4');
v.FrameRate = 10; % adjust playback speed
open(v);
%%

% 1) soft aperture (build once)
r_edge = D1/2;
rim = 0.05 * r_edge;   % 2-10% typical
Rgrid = sqrt(X.^2 + Y.^2);
cmask_soft = zeros(size(Rgrid));
inside = Rgrid <= (r_edge - rim);
transition = (Rgrid > (r_edge - rim)) & (Rgrid < (r_edge + rim));
cmask_soft(inside) = 1;
t = (Rgrid(transition) - (r_edge - rim)) / (2*rim);
cmask_soft(transition) = 0.5*(1 + cos(pi * t));
% outside remains 0

num_round_trips = 100;
%E = E.*cmask1; % clip the beam before it leaves
P0 = sum(sum(abs(E).^2)); % initial power

for i = 1:num_round_trips
    [step, Z_traveled, Z_position, E, Es] = R_L(step, Z_traveled, Z_position, E, Es, save_interval, num_steps, dz, L, H, R, amask, cmask_soft);
    E = E.*cmask2.*rmask2.*tmask;
    [step, Z_traveled, Z_position, E, Es] = L_R(step, Z_traveled, Z_position, E, Es, save_interval, num_steps, dz, L, H, R, amask, cmask_soft);
    E = E.*cmask1.*rmask1.*tmask;

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
    set(gcf,'Color','w'); % white figure background
    set(gca,'Color','w'); % white axes background
    axis('square')
    axis tight;
    view(2) % 2D view
    getframe();

    % Capture frame
    frame = getframe(gcf);
    writeVideo(v, frame);

    % Rescale 
    Pk(i) = max(max(abs(E).^2)); % peak intensity at this point
    E = E*sqrt(P0/sum(sum(abs(E).^2))); % rescale E
    BeamWidth(i) = sum(sum(sqrt(X.^2+Y.^2).*abs(E).^2))/P0;
end

close(v); % save video

