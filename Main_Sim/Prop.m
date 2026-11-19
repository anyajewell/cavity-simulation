function [laser, outputs] = Prop(consts, sim, laser, frame, outputs, toggles, dz_step)

    z2 = laser.pos + dz_step;
    z1 = laser.pos;

    % Propagation using angular spectrum method (Schmidt et al.)
    [~, laser.Gau] = ang_spec_prop(laser.Gau,laser.Ld,sim.dx,sim.dx,abs(dz_step));

    if frame.accel ~= 0
        dt_step = abs(dz_step)/consts.c;
        dv = frame.accel*dt_step; % step in velocity
        dTheta = dv/consts.c; % tilt due to transverse acceleration
        laser.Gau = laser.Gau .* exp(-1i * laser.k0 * dTheta.*sim.X); % transverse acceleration tilt phase screen
    end
    
    if frame.Omega ~= 0
        % New X-position after trajectory along characteristic curve
        rot_shift = 0.5*frame.Omega/consts.c * (z2.^2-z1.^2);
      
        % Analytic integration of phase term along characteristic 
        Xint = sim.X*(z2-z1) - 0.5*frame.Omega/consts.c * ...
            (1/3 * (z2.^3-z1.^3) + 0.5*z1*(z2-z1));

        laser.Gau = laser.Gau .* exp(-1i * laser.k0 * frame.Omega / consts.c * Xint); % rotational tilt phase screen
        
        if rot_shift ~= 0

            % Implementation of shift interpolation
            %x = sim.X(1,:).';
            %laser.Gau = interp1(x, laser.Gau.', x + rot_shift, 'spline', 'extrap').';
            % Old way:
            Xnew = sim.X - rot_shift;
            laser.Gau = interp2(Xnew,sim.Y,laser.Gau,sim.X,sim.Y,'spline'); % interpolate onto the new grid

        end

    end

    if toggles.absorbing_mask == true
        laser.Gau = laser.Gau.*sim.mask_abs;
    end
    
    if toggles.track_centers == true && toggles.outputs_switch == true
        outputs.centerx(end+1) = trapz(trapz(sim.X.*abs(laser.Gau).^2))/trapz(trapz(abs(laser.Gau).^2)); % track center x
        outputs.centery(end+1) = trapz(trapz(sim.Y.*abs(laser.Gau).^2))/trapz(trapz(abs(laser.Gau).^2)); % track center y
    end

    if toggles.outputs_switch == true && strcmp(toggles.videoplot_frequency, 'every step')
        % Plot every step here
        outputs = Write_Video_Frame(sim, laser, toggles, outputs);
    end

    if toggles.outputs_switch == true
        outputs.zs(end+1) = laser.pos;
    end

end

