function [consts, sim, laser, frame, mirror, outputs, toggles] = Initialize_Sim()
    % Constants
    c = 3e8; % [m/s]
    eps0 = (1/(36*pi))*10^(-9); % vacuum permittivity, [F/m]
    consts.c = c; consts.eps0 = eps0;
    
    % Grid
    Lx = 4; % Length of square transverse domain (one side), [m]
    N = 511; % sampling number
    dx = Lx/N; % step size 
    sim.Lx = Lx; sim.N = N; sim.dx = dx;
    
    % Initial conditions
    Ld = 1064*1e-9; % Laser wavelength, [m]
    w0 = 5e-1; % initial beam width, [m]
    k0 = (2*pi)/Ld; % wavenumber
    laser.Ld = Ld; laser.w0 = w0; laser.k0 = k0;
    
    % Frame settings
    Omega = 0.001; % rotational velocity, [rad/sec]
    accel = 1e-5; % transverse acceleration, [m/s^2]
    v0 = 0; % starting transverse velocity, [m/s]
    frame.Omega = Omega; frame.accel = accel; frame.v0 = v0;
    
    % Simulation settings
    Zmax = 75e3; % destination of z, [m] 
    Z0 = -Zmax; % starting location, [m]
    L = abs(Zmax - Z0); % cavity length, [m]
    Nz = 1e2; % number of steps in one pass across the cavity (1/2 a round trip)
    dz = L/Nz; % step size, [m]
    t0 = 0; t(1) = t0;
    dt = dz/c; % time step, [s]
    tmax = dt*Nz; % max time, [s]
    RTs = 10; % number of round trips to take
    sim.Zmax = Zmax; sim.Z0 = Z0; sim.L = L; sim.Nz = Nz; sim.dz = dz; sim.t0 = t0; sim.t = t; sim.dt = dt; sim.tmax = tmax; sim.RTs = RTs;
    
    % Mirrors
    D1 = 1; % large size to reduce clipping, [m]
    D2 = D1; % diameter of mirror 2, [m]
    Rc1 = L; % radius of curvature for mirror 1, [m]
    Rc2 = Rc1; % radius of curvature for mirror 2, [m]
    mirror(1).D = D1; mirror(2).D = D2; mirror(1).Rc = Rc1; mirror(2).Rc = Rc2;
    
    % Calculations
    zr = pi*w0^2/Ld; % Rayleigh range
    wz = w0*sqrt(1+(L/zr).^2); % spot size at mirror, analytic solution
    
    % Set up video
    videoname = sprintf('Omega=%.3f_accel=%.2e_L=%.0f_D=%.2f_RTs=%.0f.mp4', Omega, accel, L, D1, RTs);
    [saveFolder, v] = Set_Up_Video(videoname);
    outputs.saveFolder = saveFolder; outputs.v = v;
    
    % Set up simulation
    [x, y, X, Y] = Create_Grid(N, Lx, dx);
    sim.x = x; sim.y = y; sim.X = X; sim.Y = Y;
    Gau_ini = (1/(w0*pi*0.5))*exp(-(X.^2+Y.^2)./(w0^2));
    z = linspace(Z0,Zmax,Nz); % z locations to evaluate field
    z2 = linspace(Zmax,Z0,Nz); % other half, going other way
    RT = [z, z2]; % represents all positions being evaluated within one round-trip
    zs = [];
    for i = 1:RTs
        zs = [zs, RT];
    end
    sim.z = z; sim.zs = zs;
    t0 = 0;
    t = linspace(t0,tmax,Nz); % timesteps
    Gau = Gau_ini;
    centerx(1) = trapz(trapz(X.*abs(Gau).^2))/trapz(trapz(abs(Gau).^2));
    centery(1) = trapz(trapz(Y.*abs(Gau).^2))/trapz(trapz(abs(Gau).^2));
    loss_frac = zeros(1, RTs); % pre-allocate, to be filled
    laser.Gau_ini = Gau_ini; laser.Gau = Gau; outputs.centerx = centerx; outputs.centery = centery; outputs.loss_frac = loss_frac;
     
    % Mirror masks: reflecting lens phase screens and clipping masks
    rmask1 = exp(1i*k0*(X.^2+Y.^2)/(Rc1)); % reflection mask mirror 1 (RHS)
    rmask2 = exp(1i*k0*(X.^2+Y.^2)/(Rc2)); % reflection mask mirror 2 (LHS)
    cmask1 = (X.^2 + Y.^2 <= (D1/2)^2); % clipping mask mirror 1 (RHS)
    cmask2 = (X.^2 + Y.^2 <= (D2/2)^2); % clipping mask mirror 2 (LHS)
    mirror(1).rmask = rmask1; mirror(2).rmask = rmask2; mirror(1).cmask = cmask1; mirror(2).cmask = cmask2;
    
    % Settings 
    track_centers = true;
    toggles.track_centers = true;
end