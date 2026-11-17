function [EminusNext, EplusAtLHS, outputs] = ...
    L_3D(consts, sim, laser, frame, mirror, outputs, toggles)

    Nz = sim.Nz;
    mirror(2) = mirror(2);
    Eminus = laser.Eminus;
    EminusNext = zeros(size(Eminus), 'like', Eminus);
    
    slice = struct('Ld', laser.Ld, 'k0', laser.k0, ...
        'Gau', [], 'pos', []);
    
    % Record diagnostics at the mirror plane.
    propToggles = toggles;
    propToggles.outputs_switch = false;
    
    % Advance interior slices one cell to the left.
    for j = 2:Nz
        slice.Gau = Eminus(:,:,j);
        slice.pos = sim.z(j);
        slice = Prop(consts, sim, slice, frame, struct(), ...
            propToggles, -sim.dz);
        EminusNext(:,:,j-1) = slice.Gau;
    end
    
    % Travel from the first cell center to the LHS mirror.
    incidentToggles = propToggles;
    
    slice.Gau = Eminus(:,:,1);
    slice.pos = sim.z(1);
    dzToMirror = mirror(2).loc - slice.pos;
    slice = Prop(consts, sim, slice, frame, struct(), ...
        incidentToggles, dzToMirror);
    E_incident = slice.Gau;
    
    % Interact with mirror 2 (LHS)
    phi_tilt = 2*laser.k0*(mirror(2).dtheta_x*sim.X + mirror(2).dtheta_y*sim.Y); % query mirror misalignment
    tilt_mask = exp(1i*phi_tilt);
    E_reflected = E_incident .* mirror(2).cmask .* mirror(2).rmask .* tilt_mask; % clip and shape the beam
    
    if toggles.outputs_switch
        pixel_area = sim.dx*sim.dy;
        RP_before = pixel_area*sum(abs(E_incident).^2, 'all'); % power before mirror 2
        RP_after = pixel_area*sum(abs(E_reflected).^2, 'all'); % power after mirror 2
    
        reflectance = NaN;
        if RP_before > 0
            reflectance = RP_after/RP_before; % incident over reflected power
        end
    
        outputs.R2(end+1) = reflectance;
        outputs.loss2(end+1) = 1 - reflectance;
        outputs.Imax(end+1) = 0.5*consts.c*consts.eps0 ...
            * max(abs(E_reflected).^2, [], 'all');
    end
    
    % Return to the first cell center as a right-going field.
    slice.Gau = E_reflected;
    slice.pos = mirror(2).loc;
    dzFromMirror = sim.z(1) - slice.pos;
    slice = Prop(consts, sim, slice, frame, struct(), ...
        propToggles, dzFromMirror);
    
    EplusAtLHS = slice.Gau;

end