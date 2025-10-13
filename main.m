clc; clear all;

% --- SETUP --- %

% Constants
consts.c = 3*10^8; % speed of light, [m/s]
consts.eps0 = (1/(36*pi))*10^(-9); % vacuum permittivity, [F/m]

% Adjustable parameters
L = 400; % length of cavity, [m]
%D1 = 0.0254/2; % diameter of mirror 1, [m]
D1 = 0.0254;
D2 = D1; % diameter of mirror 2, [m]
Rc1 = 800; % radius of curvature for mirror 1, [m]
Rc2 = Rc1; % radius of curvature for mirror 2, [m]
N = 2048; % number of mesh points along each dim of mesh grid
lambda = 1.064e-6; % laser wavelength, [m]
W = 8*D1; % domain half width, [m]
CFL = 0.0625; % CFL number
Omega = 1; % relative rotation of spacecraft frame to inertial geocentric frame, [rad/s]

% Grid
k0 = 2*pi/lambda; % freespace wavenumber, [m^-1]
x = linspace(-W,W,N);
y = x;
dx = x(2) - x(1);
dy = dx;
dz = CFL*4*k0*dx^2; % CFL-like condition, [m]
dz = L; % make each step a trip across the cavity, [m]
dz = 0.1;
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
% --- Absorbing mask B: radial profile smoothed with separable Gaussian convolution ---
% Tunable params
taper_width_m = 0.30 * W;   % make this large (0.25-0.40 * W)
inner_radius = W - taper_width_m;
outer_radius = W;           % full absorb at W
poly_power = 2;             % polynomial fall (1 linear, 2 quadratic, etc.)

% Make base mask (polynomial roll)
ra = sqrt(X.^2 + Y.^2);
amask_base = ones(size(ra));
zone = (ra > inner_radius) & (ra < outer_radius);
s = (ra(zone) - inner_radius) ./ (outer_radius - inner_radius); % 0..1
amask_base(zone) = (1 - s).^poly_power;
amask_base(ra >= outer_radius) = 0;
amask = 1;

% Gaussian smoothing by 1D separable kernel (cheap & robust)
% choose sigma in grid points (not meters). More sigma => smoother.
sigma_m = 0.08 * W;                      % smoothing width in meters (try 0.05-0.12*W)
sigma_px = max(1, round(sigma_m / dx));  % convert to pixels
% Create 1D Gaussian kernel
k_half = ceil(4 * sigma_px);
xk = -k_half:k_half;
g1d = exp( - (xk.^2) / (2 * sigma_px^2) );
g1d = g1d / sum(g1d); % normalize

% separable convolution: first along x, then y
amask_sm = conv2(g1d, g1d, amask_base, 'same');

% enforce bounds and numerical floor
amask_sm(amask_sm<1e-12) = 0;
amask = amask_sm;

% visualize
figure(101); imagesc(x, y, amask); axis equal tight; colorbar;
title(sprintf('Smoothed mask (taper %.3fm, sigma_px=%d)', taper_width_m, sigma_px));
drawnow;

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

fileName = sprintf('Omega=%.3f_L=%.0fkm_Rc=%.0f.mp4', Omega, L, Rc1);
filePath = fullfile(saveFolder, fileName);
v = VideoWriter(filePath, 'MPEG-4');
v.FrameRate = 10; % adjust playback speed
open(v);

%% SIMULATION, R --> L, one trip across in small steps
% Just go once across the cavity

P0 = sum(sum(abs(E).^2)); % initial power

% Clip and shape the beam
E = E.*cmask1.*rmask1.*tmask; % beam leaves from the RHS (mirror1)

for n = 1:num_steps
    step = step + 1;
    Z_traveled(step) = step*dz; % record total distance propagated so far

    % Record current position within the cavity, AFTER the step has
    % been taken
    if step == 1 % this is the first step
        Z_position(step) = L/2 - dz;
    else
        Z_position(step) = Z_position(step-1) - dz;
    end

    % Propagation
    FE = fft2(E); % transform beam to frequency domain
    FE = FE.*fftshift(H).*R(Z_position(step), Z_position(step)-dz); % propagate beam in frequency domain 
    E = ifft2(FE); % transform back to space domain 
    E = tmask.*E.*amask; % apply the tilting mask
    
    % Save E field snapshots and write video frame
    if mod(step, save_interval) == 0
        step_label = sprintf('step_%d', step);
        Es.(step_label) = E; % save intermediate field

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
    ylim([-2*D1, 2*D1])
    xlim([-2*D1, 2*D1])
    view(2) % 2D view
    getframe();

    % Capture frame
    frame = getframe(gcf);
    writeVideo(v, frame);
    end

    % Rescale 
    % Pk(n) = max(max(abs(E).^2)); % peak intensity at this point
    % E = E*sqrt(P0/sum(sum(abs(E).^2))); % rescale E
    % BeamWidth(n) = sum(sum(sqrt(X.^2+Y.^2).*abs(E).^2))/P0;
end

E = E.*cmask2.*rmask2.*tmask; % beam arrives at the LHS (mirror2)
close(v); % save video

%% SIMULATION

% Set how many round trips across the cavity you want to take.
% e.g. RHS --> LHS --> RHS = 1 round trip.

num_round_trips = 100;
E = E.*cmask1; % clip the beam before it leaves
P0 = sum(sum(abs(E).^2)); % initial power

for i = 1:num_round_trips
    [step, Z_traveled, Z_position, E, Es] = R_L(step, Z_traveled, Z_position, E, Es, save_interval, num_steps, dz, L, H, R, amask, consts);
    E = E.*cmask2.*rmask2.*tmask;
    [step, Z_traveled, Z_position, E, Es] = L_R(step, Z_traveled, Z_position, E, Es, save_interval, num_steps, dz, L, H, R, amask, consts);
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
    ylim([-2*D1, 2*D1])
    xlim([-2*D1, 2*D1])
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

%% Post-processing

% Trackers used for both numeric and analytic
dz_post = 1; % little steps across the cavity  
zs = -L/2 : dz_post : +L/2; % z positions, measured relative to the cavity center, z = 0, [m]
numSteps = length(zs); % number of steps to be taken for one trip across the cavity

% Numeric
H_post = exp(1i/(2*k0)*dz_post*(KX.^2+KY.^2)); % compute free space transfer function for a shorter distance
E0_num = E; % Use the last saved E, numeric complex field at current plane   
Es_num = cell(1, numSteps); wzs = zeros(1, numSteps); % initialize cell arrays to save numeric results
[Es_num, wzs] = Save_Numeric_Results(E0_num, Es_num, wzs, numSteps, H_post, X, Y, dx, dy);  
waist_index = find(wzs == min(wzs), 1); % select where the numeric waist occurs
I = abs(Es_num{waist_index}).^2; % numeric intensity at waist
I0 = max(I(:)); % on-axis peak intensity
R = sqrt(X.^2 + Y.^2); % radial coordinate grid
[theta, rvals] = cart2pol(X, Y); % radially average intensity to smooth noise/asymmetry
dr = dx; rvec = 0:dr:max(R(:)); % incremental radii of the beam, for options
Iradial = zeros(size(rvec));

for k = 1:length(rvec)-1
    mask = (R >= rvec(k)) & (R < rvec(k+1));
    Iradial(k) = mean(I(mask)); % symmetric intensity profile
end

target = I0*exp(-2); % target is 1/e^2 of the max instensity, can set different target
[~,idx] = min(abs(Iradial - target));
w0_num = rvec(idx); % numeric beam waist
fprintf('Numeric w0: %.3g mm \n', w0_num*1e3);
zr_num = pi * w0_num^2 / lambda; % numeric Rayleigh range
wz_num = w0_num * sqrt(1+(zs./zr_num).^2); % numeric spot size w(z)


% Analytic 
R = Rc1; % set which mirror the beam is arriving at
w0_ana = sqrt(lambda*L/pi * sqrt(g1*g2*(1-g1*g2)/(g1+g2-2*g1*g2)^2)); % analytic beam waist from cavity geometry
fprintf('Analytic w0: %.3g mm \n', w0_ana*10^3);
zr = pi * w0_ana^2 / lambda; % Rayleigh range from the analytic waist
wz_ana = w0_ana * sqrt(1+(zs./zr).^2); % analytic spot size along zs
q0 = 1i*(pi*w0_ana^2)/lambda; % initial complex radius, evaluated at the center of the cavity, i.e. at the beam waist
qz = R - (pi*wz_ana.^2)/(1i*lambda); % complex radius of curvature
%E_ana = 1/q0 * exp(-1i*k0 * (X.^2+Y.^2)/(2*q0)); % analytical field at beam waist

%% Graphing

figure;
plot(zs, +wz_ana, 'b', 'LineWidth', 1.5, 'DisplayName', 'Analytic w(z)'); hold on; 
plot(zs, -wz_ana, 'b', 'LineWidth', 1.5, 'HandleVisibility', 'off');
plot(zs, +wz_num, 'r', 'LineWidth', 1.5, 'DisplayName', 'Numeric w(z)'); hold on; 
plot(zs, -wz_num, 'r', 'LineWidth', 1.5, 'HandleVisibility', 'off');
xlabel('z [m]');
ylabel('Transverse position [m]');
title('Beam cross section');
grid on;
legend()
axis tight;


%% Graveyard
