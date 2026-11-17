function [laser, outputs, sim, gain_medium] = R(consts, sim, laser, frame, mirror, outputs, toggles, gain_medium)        
    
    while laser.pos < mirror(1).loc
        frame.Omega = frame.Omega; % opportunity to read/update Omega from other code
        frame.accel = frame.accel; % opportunity to read/update accel from other code
        if laser.pos == 0 % beam at the center of the cavity
            outputs.Gau_center = laser.Gau; % save center mode
        end
        dt_to_mirror = (mirror(1).loc - laser.pos)/consts.c;
        dt_step = min(sim.dt, dt_to_mirror);
        [laser, outputs] = Prop(consts, sim, laser, frame, outputs, toggles, dt_step);
        laser.pos = laser.pos + consts.c*dt_step;
        if dt_step == dt_to_mirror
            laser.pos = mirror(1).loc;
        end
    end

    RP_before = sim.dx^2*sum(abs(laser.Gau).^2,'all'); % power before mirror 1

    % Set up tilt mask
    theta_x1 = mirror(1).dtheta_x; theta_y1 = mirror(1).dtheta_y; % query misalignment
    phi_tilt1 = 2*laser.k0*(theta_x1*sim.X + theta_y1*sim.Y);
    mirror(1).tmask = exp(1i*phi_tilt1); % tilt mask

    % Interact with mirror 1
    laser.Gau = laser.Gau.*mirror(1).cmask.*mirror(1).rmask.*mirror(1).tmask; % clip and shape the beam
    
    RP_after = sim.dx^2*sum(abs(laser.Gau).^2,'all'); % power after mirror 1
    I_after = 0.5*consts.c*consts.eps0*abs(laser.Gau).^2; % intensity after mirror 1
    laser.P_ref(end+1) = RP_after; % reference power for RT loss calculation 
    
    if toggles.outputs_switch == true
        outputs.R1(end+1) = RP_after / RP_before; % reflected over incident power
        outputs.loss1(end+1) =  1 - RP_after / RP_before;
        outputs.Imax(end+1) = max(I_after,[],'all');
    end

    % Apply gain
    if toggles.gain_switch == true % Gain medium present at mirror 1
        P_before = sim.dx*sim.dx*sum(abs(laser.Gau).^2,'all');
        I = 0.5*consts.c*consts.eps0*abs(laser.Gau).^2; % laser intensity profile
        gain_medium.g = gain_medium.g0_profile ./ (1 + I./gain_medium.I_sat); % gain function
        laser.Gau = laser.Gau.*exp(gain_medium.g/2); % rescale electric field
        P_after = sim.dx*sim.dx*sum(abs(laser.Gau).^2,'all');
        if toggles.outputs_switch == true
            outputs.gain(end+1) = P_after / P_before;
        end
    end
        
    [sim] = Turn_Around(sim);

end