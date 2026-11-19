function [EplusNext, EminusAtRHS, outputs, gain_medium] = R_3D(consts, sim, laser, frame, mirror, outputs, toggles, gain_medium)        

    Nz = sim.Nz;
    z = sim.z;
    dz = sim.dz;
    Eplus = laser.Eplus;
    EplusNext = zeros(size(Eplus), 'like', Eplus);
    
    slice = struct('Ld', laser.Ld, 'k0', laser.k0, ...
        'Gau', [], 'pos', []);
    slice_template = slice;
    
    % Record one set of mirror diagnostics per time step:
    propToggles = toggles;
    propToggles.outputs_switch = false;
    
    % Advance each interior slice one cell to the right:
    for j = 1:Nz-1
        local_slice = slice_template;
        local_slice.Gau = Eplus(:,:,j); % picks out the longitudinal cell
        local_slice.pos = z(j); % corresponding position
        local_slice = Prop(consts, sim, local_slice, frame, struct(), propToggles, dz); % propagates the slice
        EplusNext(:,:,j+1) = local_slice.Gau; % fills out the updated array
    end
    
    % Advance the last right-going slice to the RHS mirror.
    incidentToggles = propToggles;
    
    slice.Gau = Eplus(:,:,Nz); % the final slice, at the RHS of the cavity
    slice.pos = sim.z(Nz);
    dzToMirror = mirror(1).loc - slice.pos; % half step, +dz/2
    slice = Prop(consts, sim, slice, frame, struct(), ...
        incidentToggles, dzToMirror);
    E_incident = slice.Gau;
    
    % Interact final slice with mirror 1
    phi_tilt = 2*laser.k0*(mirror(1).dtheta_x*sim.X + mirror(1).dtheta_y*sim.Y); % query misalignment and apply tilt
    tilt_mask = exp(1i*phi_tilt);
    E_reflected = E_incident .* mirror(1).cmask .* mirror(1).rmask .* tilt_mask; % clip and shape the beam
    
    % Record mirror loss and intensity before applying gain.
    pixel_area = sim.dx*sim.dy;
    intensity_factor = 0.5*consts.c*consts.eps0;
    power_factor = intensity_factor*pixel_area;
    
    if toggles.outputs_switch
        RP_before = power_factor*sum(abs(E_incident).^2, 'all');
        RP_after = power_factor*sum(abs(E_reflected).^2, 'all');
    
        reflectance = NaN;
        if RP_before > 0
            reflectance = RP_after/RP_before; % reflected over incident power
        end
    
        outputs.R1(end+1) = reflectance;
        outputs.loss1(end+1) = 1 - reflectance;
        outputs.Imax(end+1) = intensity_factor ...
            * max(abs(E_reflected).^2, [], 'all');
    end
    
    % Apply the existing instantaneous gain screen after reflection.
    if toggles.gain_switch % Gain medium present at mirror 1
        I = intensity_factor*abs(E_reflected).^2;
        gain_medium.g = gain_medium.g0_profile ./ ...
            (1 + I./gain_medium.I_sat); % gain function
        E_reflected = E_reflected .* exp(gain_medium.g/2); % rescale electric field
    
        if toggles.outputs_switch
            gain_ratio = NaN;
            if RP_after > 0
                RP_after_gain = power_factor*sum(abs(E_reflected).^2, 'all');
                gain_ratio = RP_after_gain/RP_after; % gain
            end
            outputs.gain(end+1) = gain_ratio;
        end
    end

    if toggles.outputs_switch
        outputs.Gau_RHS_mirror = E_reflected; % save for comparison with 2D solver
    end
    
    % Outermost, reflected slice becomes a left-going field
    slice.Gau = E_reflected;
    slice.pos = mirror(1).loc;
    dzFromMirror = sim.z(Nz) - slice.pos; % half step, -dz/2
    slice = Prop(consts, sim, slice, frame, struct(), ...
        propToggles, dzFromMirror);
    
    EminusAtRHS = slice.Gau;

end