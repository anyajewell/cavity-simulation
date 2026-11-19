function laser = Fill_Cavity(consts, sim, laser, frame, mirror, toggles, gain_medium)

    Eseed = laser.Gau;

    % Filled cavity initial condition
    % Start with the right-going field at the first cell center.
    [Ny, Nx] = size(Eseed);
    laser.Eplus  = zeros(Ny, Nx, sim.Nz, 'like', Eseed);
    laser.Eminus = zeros(Ny, Nx, sim.Nz, 'like', Eseed);
    laser.Eplus(:,:,1) = Eseed;
    
    % Build the initial field without recording propagation diagnostics.
    initToggles = toggles;
    initToggles.outputs_switch = false;
    
    slice = struct('Ld', laser.Ld, 'k0', laser.k0, ...
        'Gau', Eseed, 'pos', sim.z(1));
    
    % Fill the right-going field from left to right.
    for j = 1:sim.Nz-1
        slice.pos = sim.z(j);
        dz_step = sim.z(j+1) - sim.z(j);
        slice = Prop(consts, sim, slice, frame, struct(), initToggles, dz_step);
        laser.Eplus(:,:,j+1) = slice.Gau;
    end
    
    % Use the existing RHS reflection to seed the left-going field
    % at the rightmost cell center.
    [~, EminusAtRHS] = R_3D(consts, sim, laser, frame, mirror, ...
        struct(), initToggles, gain_medium);
    
    laser.Eminus(:,:,sim.Nz) = EminusAtRHS;
    slice.Gau = EminusAtRHS;
    
    % Fill the left-going field from right to left.
    for j = sim.Nz:-1:2
        slice.pos = sim.z(j);
        dz_step = sim.z(j-1) - sim.z(j);
        slice = Prop(consts, sim, slice, frame, struct(), initToggles, dz_step);
        laser.Eminus(:,:,j-1) = slice.Gau;
    end

end