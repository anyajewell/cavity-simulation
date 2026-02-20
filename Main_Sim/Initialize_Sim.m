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
    
    % Interlink frame settings
    Omega = 0.001; % rotational velocity, [rad/sec]
    accel = 1e-5; % transverse acceleration, [m/s^2]
    v0 = 0; % starting transverse velocity, [m/s]
    frame.Omega = Omega; frame.accel = accel; frame.v0 = v0;
    
    % Simulation settings
    %Zmax = 75e3; % destination of z, [m] 
    Z0 = 0; % starting location, arbitrary, anywhere within the cavity [m]
    L = 150e3; % cavity length, [m]
    Nz = 1e2; % number of steps in one pass across the cavity (1/2 a round trip)
    dz = L/Nz; % step size, sign determines initial direction [m]
    t0 = 0; t(1) = t0;
    dt = dz/c; % time step, [s]
    tmax = dt*Nz; % max time, [s]
    RTs = 50; % number of round trips to take
    sim.Zmax = Zmax; sim.Z0 = Z0; sim.L = L; sim.Nz = Nz; sim.dz = dz; sim.t0 = t0; sim.t = t; sim.dt = dt; sim.tmax = tmax; sim.RTs = RTs;
    
    % Mirrors
    D1 = 1; % large size to reduce clipping, [m]
    D2 = D1; % diameter of mirror 2, [m]
    Rc1 = L; % radius of curvature for mirror 1, [m]
    Rc2 = Rc1; % radius of curvature for mirror 2, [m]
    loc1 = L/2; % location of mirror 1, RHS of cavity [m]
    loc2 = -L/2; % location of mirror 2, LHS of cavity [m]
    mirror(1).D = D1; mirror(2).D = D2; mirror(1).Rc = Rc1; mirror(2).Rc = Rc2; mirror(1).loc = loc1; mirror(2).loc = loc2;
    
    % Calculations
    zr = pi*w0^2/Ld; % Rayleigh range
    wz = w0*sqrt(1+(L/zr).^2); % spot size at mirror, analytic solution
    theta_D = Ld/D1; % diffraction angle
    N_F = (D1/2)^2 / (L*Ld); % Fresnel number (>>1 means diffraction losses are small and higher order modes can be supported)
    laser.zr = zr; laser.wz = wz; laser.theta_D = theta_D; laser.N_F = N_F; laser.pos = Z0;
    
    % Misalignment angles
    theta_x1 = 0.05*theta_D; % misalignment based on diffraction angle, [rad]
    sim.theta_x1 = theta_x1;

    % Set up video
    videoname = sprintf('Omega=%.3f_accel=%.2e_L=%.0f_D=%.2f_RTs=%.0f.mp4', Omega, accel, L, D1, RTs);
    [saveFolder, v] = Set_Up_Video(videoname);
    outputs.saveFolder = saveFolder; outputs.v = v;
    
    % Set up simulation
    [x, y, X, Y] = Create_Grid(N, Lx, dx);
    sim.x = x; sim.y = y; sim.X = X; sim.Y = Y;
    Gau_ini = (1/(w0*pi*0.5))*exp(-(X.^2+Y.^2)./(w0^2));
    Pseed = 1; % choose laser seed power, [W]
    Gau_ini = Normalize_Laser_Field_To_Power(Gau_ini, Pseed, sim.dx, sim.dx, consts.c, consts.eps0); % scale profile to laser power
    
    % Set up zs
    z_LR = linspace(mirror(2).loc, mirror(1).loc, Nz); % z locations to evaluate field from LHS to RHS
    z_RL = linspace(mirror(1).loc, mirror(2).loc, Nz); % other half, going other way from RHS to LHS

    % Determine how to structure 
    if Z0 ~= mirror(1).loc && Z0 ~= mirror(2).loc % wavefront is not starting at either end of the cavity
        if sign(dz) > 0 % laser traveling right
            z1 = 
        else % laser traveling left
            z1 = 
        end
    else
        if sign(dz) > 0 % laser traveling L --> R
            z1 = z_LR;
        else % laser traveling R --> L
            z1 = z_RL;
        end
    end
    
    RT = [z1, z2]; % represents all positions being evaluated within first trip
    zs = []; zs(1) = z1;

    for i = 1:RTs
        zs = [zs, RT]; % represents ALL positions, in order, this wavefront will be evaluated at
    end

    sim.z = z_LR; sim.zs = zs;

    t0 = 0;
    t = linspace(t0,tmax,Nz); % timesteps
    Gau = Gau_ini;
    centerx(1) = trapz(trapz(X.*abs(Gau).^2))/trapz(trapz(abs(Gau).^2));
    centery(1) = trapz(trapz(Y.*abs(Gau).^2))/trapz(trapz(abs(Gau).^2));
    loss_frac = zeros(1, RTs); % pre-allocate, to be filled
    loss1 = []; loss2 = []; R1 = []; R2 = []; gain = []; power = [];
    laser.Gau_ini = Gau_ini; laser.Gau = Gau; outputs.centerx = centerx; outputs.centery = centery; outputs.loss_frac = loss_frac;
    outputs.loss1 = loss1; outputs.loss2 = loss2; outputs.R1 = R1; outputs.R2 = R2; outputs.gain = gain; outputs.Imax = power;
     
    % Mirror masks: reflecting lens phase screens and clipping masks
    rmask1 = exp(1i*k0*(X.^2+Y.^2)/(Rc1)); % reflection mask mirror 1 (RHS)
    rmask2 = exp(1i*k0*(X.^2+Y.^2)/(Rc2)); % reflection mask mirror 2 (LHS)
    cmask1 = (X.^2 + Y.^2 <= (D1/2)^2); % clipping mask mirror 1 (RHS)
    cmask2 = (X.^2 + Y.^2 <= (D2/2)^2); % clipping mask mirror 2 (LHS)
    mirror(1).rmask = rmask1; mirror(2).rmask = rmask2; mirror(1).cmask = cmask1; mirror(2).cmask = cmask2;
    
    % Settings 
    track_centers = true;
    gain_switch = false; % gain ON or OFF
    toggles.track_centers = true;
    toggles.gain_switch = gain_switch;
end