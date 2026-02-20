function [laser, outputs, sim, gain_medium] = R(consts, sim, laser, frame, mirror, outputs, toggles, gain_medium, a)        
    
    while laser.pos < mirror(1).loc
        frame.Omega = frame.Omega; % opportunity to read/update Omega from other code
        frame.accel = frame.accel; % opportunity to read/update accel from other code
        if laser.pos + sim.dz > mirror(1).loc % beam will step past mirror
            dz_step = mirror(1).loc - laser.pos; % calculate fractional step
            if dz_step == 0 % beam will step directly to mirror 1
                if toggles.track_centers == true
                    outputs.centerx(end+1) = trapz(trapz(sim.X.*abs(laser.Gau).^2))/trapz(trapz(abs(laser.Gau).^2)); % track center x
                    outputs.centery(end+1) = trapz(trapz(sim.Y.*abs(laser.Gau).^2))/trapz(trapz(abs(laser.Gau).^2)); % track center y
                end
                break
            else % propagate by fractional step to mirror 1
            [laser, outputs] = Prop(consts, sim, laser, frame, outputs, toggles, dz_step); % propagate fractional step to mirror
            laser.pos = laser.pos + dz_step; % update laser position
            end
        else % propagate by dz, as normal
            [laser, outputs] = Prop(consts, sim, laser, frame, outputs, toggles, sim.dz); % propagation loop
            laser.pos = laser.pos + sim.dz; % update laser position
        end
        outputs.zs(end+1) = laser.pos;
    end

    % Interact with mirror 1 (RHS)
    if toggles.gain_switch == true % Gain medium present at mirror 1
        P_before = sim.dx*sim.dx*sum(abs(laser.Gau).^2,'all');
        I = 0.5*consts.c*consts.eps0*abs(laser.Gau).^2; % laser intensity profile
        gain_medium.g = gain_medium.g0_profile ./ (1 + I./gain_medium.I_sat); % gain function
        laser.Gau = laser.Gau.*exp(gain_medium.g/2); % rescale electric field
        P_after = sim.dx*sim.dx*sum(abs(laser.Gau).^2,'all');
        outputs.gain(end+1) = P_after / P_before;
    end

    RP_before = trapz(trapz(abs(laser.Gau).^2)); % power before mirror 1
    I_before = 0.5*consts.c*consts.eps0*abs(laser.Gau).^2; % intensity before mirror 1
    theta_x1 = 0; theta_y1 = 0; % query mirror misalignment

    if a == 30 % introduce one time misalignment
        theta_x1 = sim.theta_x1;
    end

    phi_tilt1 = 2*laser.k0*(theta_x1*sim.X + theta_y1*sim.Y);
    mirror(1).tmask = exp(1i*phi_tilt1); % tilt mask
    laser.Gau = laser.Gau.*mirror(1).cmask.*mirror(1).rmask.*mirror(1).tmask; % clip and shape the beam
    
    RP_after = trapz(trapz(abs(laser.Gau).^2)); % power after mirror 1
    I_after = 0.5*consts.c*consts.eps0*abs(laser.Gau).^2; % intensity after mirror 1
    
    outputs.R1(end+1) = RP_after / RP_before; % reflected over incident power
    outputs.loss1(end+1) =  1 - RP_after / RP_before;
    outputs.Imax(end+1) = max(I_after,[],'all');
    
    [sim] = Turn_Around(sim);

end