function [consts, sim, laser, frame, mirror, outputs, toggles, gain_medium] = Initialize_Sim(N_F)
    
    % Settings 
    track_centers = true;
    gain_switch = false; % gain ON or OFF
    outputs_switch = true;
    videoplot_frequency = 'every mirror'; % 'every step', 'every mirror', or 'never/none'
    finish_line = 'convergence'; % 'convergence' or 'RTs'
    absorbing_mask = true;
    resize_grid = false;
    convergence_def = 'general'; % 'TEM00' or 'general'
    initial_profile = 'cavity mode'; % 'cavity mode' or 'general'
    toggles.track_centers = track_centers; toggles.gain_switch = gain_switch; toggles.outputs_switch = outputs_switch; 
    toggles.videoplot_frequency = videoplot_frequency; toggles.finish_line = finish_line; toggles.absorbing_mask = absorbing_mask;
    toggles.resize_grid = resize_grid; toggles.convergence_def = convergence_def; toggles.initial_profile = initial_profile;

    % Constants
    c = 3e8; % [m/s]
    eps0 = (1/(36*pi))*10^(-9); % vacuum permittivity, [F/m]
    consts.c = c; consts.eps0 = eps0;
    
    % Grid
    Lx = 4; % Length of square transverse domain (one side), [m]
    N = 512; % sampling number
    dx = Lx/N; % step size 
    sim.Lx = Lx; sim.N = N; sim.dx = dx;
    
    % Interlink frame settings
    Omega = 0.001; % rotational velocity, [rad/sec]
    accel = 1e-5; % transverse acceleration, [m/s^2]
    v0 = 0; % starting transverse velocity, [m/s]
    frame.Omega = Omega; frame.accel = accel; frame.v0 = v0;
    
    % Simulation settings
    L = 150e3; % cavity length, [m]
    %D1 = 1;
    Z0 = -L/2; % starting location, arbitrary, anywhere within the cavity [m]
    Nz = 1e2; % number of steps in one pass across the cavity (1/2 a round trip)
    dz = L/Nz; % step size, sign determines initial direction [m]
    t0 = 0; t(1) = t0;
    dt = abs(dz/c); % time step, [s]
    tmax = dt*Nz; % max time, [s]
    z = linspace(Z0, L, Nz); % evaluatory locations within cavity
    sim.Z0 = Z0; sim.L = L; sim.Nz = Nz; sim.dz = dz; sim.t0 = t0; sim.t = t; sim.dt = dt; sim.tmax = tmax; sim.z = z;
    RTs = 500; % number of round trips to take, user set (represents a max if finish_line 'convergence' is on)
    %RTs = Set_Max_RTs(sim, consts, sampling_time); % to be used with PCAC
    sim.RTs = RTs;
    
    % Mirrors
    Ld = 1064*1e-9; % Laser wavelength, [m]
    D1 = 2*sqrt(N_F*L*Ld);
    %D1 = .8; % large size to reduce clipping, [m]
    D2 = D1; % diameter of mirror 2, [m]
    Rc1 = L; % radius of curvature for mirror 1, [m]
    Rc2 = Rc1; % radius of curvature for mirror 2, [m]
    loc1 = L/2; % location of mirror 1, RHS of cavity [m]
    loc2 = -L/2; % location of mirror 2, LHS of cavity [m]
    dtheta_x = 0; dtheta_y = 0; % intial mirror misalignment, [rad]
    mirror(1).D = D1; mirror(2).D = D2; mirror(1).Rc = Rc1; mirror(2).Rc = Rc2; mirror(1).loc = loc1; mirror(2).loc = loc2; mirror(1).dtheta_x = dtheta_x; mirror(1).dtheta_y = dtheta_y;
    
    % Initial laser conditions
    Ld = 1064*1e-9; % Laser wavelength, [m]
    w0 = 5e-1; % initial beam width, [m]
    k0 = (2*pi)/Ld; % wavenumber
    zr = pi*w0^2/Ld; % Rayleigh range
    wz = w0*sqrt(1+(L/zr).^2); % spot size at mirror, analytic solution
    theta_D = Ld/D1; % diffraction angle
    %N_F = (D1/2)^2 / (L*Ld); % Fresnel number (>>1 means diffraction losses are small and higher order modes can be supported)
    laser.Ld = Ld; laser.w0 = w0; laser.k0 = k0; laser.zr = zr; laser.wz = wz; laser.theta_D = theta_D; laser.N_F = N_F; laser.pos = Z0;

    % Set up video
    videoname = sprintf('Omega=%.3f_accel=%.2e_L=%.0f_D=%.2f_RTs=%.0f.mp4', Omega, accel, L, D1, RTs);
    if toggles.outputs_switch == true
        [saveFolder, v] = Set_Up_Video(videoname);
        outputs.saveFolder = saveFolder; outputs.v = v;
    else
        outputs = struct();
    end
    
    % Set up simulation
    [x, y, X, Y] = Create_Grid(N, Lx, dx);
    sim.x = x; sim.y = y; sim.X = X; sim.Y = Y;

    if strcmp(initial_profile, 'general') % initialize general gaussian beam
        Gau_ini = (1/(w0*pi*0.5))*exp(-(X.^2+Y.^2)./(w0^2));
    else
        Gau_ini = load("C:\Users\Anya Jewell\Documents\GitHub\cavity-simulation\Results\Loss_vs_NF\NF=1\workspace_2026_03_17_144605.mat").outputs.Gau_LHS;
    end

    P_ref = []; % to store reference power for round-trip loss calculation at RHS
    %Pseed = 1; % choose laser seed power, [W]
    %Gau_ini = Normalize_Laser_Field_To_Power(Gau_ini, Pseed, sim.dx, sim.dx, consts.c, consts.eps0); % scale profile to laser power
    t0 = 0;
    t = linspace(t0,tmax,Nz); % timesteps
    Gau = Gau_ini;
    centerx(1) = trapz(trapz(X.*abs(Gau).^2))/trapz(trapz(abs(Gau).^2));
    centery(1) = trapz(trapz(Y.*abs(Gau).^2))/trapz(trapz(abs(Gau).^2));
    loss1 = []; loss2 = []; R1 = []; R2 = []; gain = []; Imax = []; zs = []; loss_frac = [];
    laser.Gau_ini = Gau_ini; laser.Gau = Gau; laser.P_ref = P_ref;
    if toggles.outputs_switch == true
        outputs.centerx = centerx; outputs.centery = centery; outputs.loss_frac = loss_frac;
        outputs.loss1 = loss1; outputs.loss2 = loss2; outputs.R1 = R1; outputs.R2 = R2; outputs.gain = gain; outputs.Imax = Imax; outputs.zs = zs;
    end
     
    % Mirror masks: reflecting lens phase screens and clipping masks
    rmask1 = exp(-1i*k0*(X.^2+Y.^2)/(Rc1)); % reflection mask mirror 1 (RHS), negative sign because propagation uses positive sign convention
    rmask2 = exp(-1i*k0*(X.^2+Y.^2)/(Rc2)); % reflection mask mirror 2 (LHS), negative sign because propagation uses positive sign convention
    cmask1 = (X.^2 + Y.^2 <= (D1/2)^2); % clipping mask mirror 1 (RHS)
    cmask2 = (X.^2 + Y.^2 <= (D2/2)^2); % clipping mask mirror 2 (LHS)
    mirror(1).rmask = rmask1; mirror(2).rmask = rmask2; mirror(1).cmask = cmask1; mirror(2).cmask = cmask2;

    % Domain masks: absorbing mask
    if toggles.absorbing_mask == true
        r = sqrt(sim.X.^2 + sim.Y.^2);
        r0 = 0.8*(sim.Lx/2); % damping starts at 80% radius
        w = 0.1*(sim.Lx/2); % damping width
        sim.mask_abs = ones(size(r));
        idx = r > r0;
        sim.mask_abs(idx) = exp(-((r(idx)-r0)/w).^8); % steep super-Gaussian
    end

    if toggles.gain_switch == true
        gain_medium = Initialize_Gain_Medium(sim, mirror);
    else
        gain_medium = struct();   
    end
        
    if toggles.resize_grid == true
        sim.grid0.N = sim.N;
        sim.grid0.dx = sim.dx;
        sim.grid0.Lx = sim.Lx;
        sim.grid0.x = sim.x;
        sim.grid0.y = sim.y;
        [sim, laser, mirror, gain_medium] = Increase_Domain_Size(sim, laser, mirror, toggles, 1024); % for the first pass only
    end

end