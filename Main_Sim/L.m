function [laser, outputs, sim, gain_medium] = L(consts, sim, laser, frame, mirror, outputs, toggles, gain_medium)

    while laser.pos > mirror(2).loc
        frame.Omega = frame.Omega; % opportunity to read/update Omega from other code
        frame.accel = frame.accel; % opportunity to read/update accel from other code
        dt_to_mirror = (laser.pos - mirror(2).loc)/consts.c;
        dt_step = min(sim.dt, dt_to_mirror);
        [laser, outputs] = Prop(consts, sim, laser, frame, outputs, toggles, dt_step);
        laser.pos = laser.pos - consts.c*dt_step;
        if dt_step == dt_to_mirror
            laser.pos = mirror(2).loc;
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
    
    if toggles.outputs_switch == true
        outputs.R2(end+1) = RP_after / RP_before; % reflected over incident power
        outputs.loss2(end+1) =  1 - RP_after / RP_before;
        outputs.Imax(end+1) = max(I_after,[],'all');
    end

    [sim] = Turn_Around(sim);

end