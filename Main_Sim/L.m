function [laser, outputs] = L()

    for i = 1:sim.Nz
        pos = sim.z(i); % current position
        frame.Omega = Omega; % opportunity to read/update Omega from other code
        frame.accel = accel; % opportunity to read/update accel from other code
        if pos + sim.dz < mirror(2).loc 
            dz1 = sim.Zmax-pos; 
            if dz1 == 0 % beam already at mirror
                if toggles.track_centers == true
                    outputs.centerx(end+1) = trapz(trapz(sim.X.*abs(laser.Gau).^2))/trapz(trapz(abs(laser.Gau).^2)); % track center x
                    outputs.centery(end+1) = trapz(trapz(sim.Y.*abs(laser.Gau).^2))/trapz(trapz(abs(laser.Gau).^2)); % track center y
                end
                break
            else % propagate by fractional step to L
            [laser, outputs] = Prop(consts, sim, laser, frame, outputs, toggles, i); % propagate fractional step to mirror
            end
        else % propagate by dz, as normal
            [laser, outputs] = Prop(consts, sim, laser, frame, outputs, toggles, i); % propagation loop
        end
    end

    % Interact with mirror 2 (LHS)
    RP_before = trapz(trapz(abs(laser.Gau).^2)); % power before mirror 2
    I_before = 0.5*consts.c*consts.eps0*abs(laser.Gau).^2; % intensity before mirror 2
    theta_x2 = 0; theta_y2 = 0; % query mirror misalignment
    phi_tilt2 = 2*laser.k0*(theta_x2*sim.X + theta_y2*sim.Y);
    mirror(2).tmask = exp(1i*phi_tilt2); % tilt mask
    laser.Gau = laser.Gau.*mirror(2).cmask.*mirror(2).rmask.*mirror(2).tmask; % clip and shape the beam
    RP_after = trapz(trapz(abs(laser.Gau).^2)); % power after mirror 2
    I_after = 0.5*consts.c*consts.eps0*abs(laser.Gau).^2; % intensity after mirror 2
    outputs.R2(end+1) = RP_after / RP_before; % reflected over incident power
    outputs.loss2(end+1) =  1 - RP_after / RP_before;
    outputs.Imax(end+1) = max(I_after,[],'all');

    [sim] = Turn_Around(sim);

end