function [laser, outputs, gain_medium, sim] = Propagate_n_RTs_3D(consts, sim, laser, frame, mirror, outputs, toggles, gain_medium)

    if toggles.resize_grid == true
        [sim, laser, mirror, gain_medium] = Decrease_To_Grid0_Centroid(sim, laser, mirror, outputs, toggles);
    end
    
    Nz = sim.Nz;
    steps_per_RT = 2*Nz; % number of spatial steps per round-trip
    power_factor =  0.5*consts.c*consts.eps0*sim.dx*sim.dy;
    state = []; % initialize for convergence checks
    check_convergence = strcmp(toggles.finish_line, 'convergence'); % if true, we are checking convergence every round-trip
    previous_loss = NaN;

    outputs.t = sim.t;
    outputs.Pplus = power_factor*reshape( ...
        sum(sum(abs(laser.Eplus).^2, 1), 2), 1, []);
    outputs.Pminus = power_factor*reshape( ...
        sum(sum(abs(laser.Eminus).^2, 1), 2), 1, []);

    % Set up and start video
    [video, frameStride, fig, layout, sliceSurfaces, ix, iy, sliceIndices, intensityFactor] = Set_Up_3D_Video(sim,consts,laser);
    Write_Cavity_Frame(video, fig, layout, sliceSurfaces, laser, ix, iy, sliceIndices, intensityFactor, sim.t);

    % Outer Loop:
    for a = 1:sim.RTs % go for the prescribed number of round-trips

        E_start = laser.Eplus(:,:,1); % field before the round-trip begins, for comparison later, to compare the same direction at the same plane one round-trip apart
        E_rhs = [];
        rhs_reflectance = NaN;

        % Inner Loop:
        for step = 1:steps_per_RT % step through this round-trip, every step represents a spatial and time step

            [EplusNext, EminusAtRHS, outputs, gain_medium] = R_3D(consts, sim, laser, frame, mirror, outputs, toggles, gain_medium); % forward-propagating field

            [EminusNext, EplusAtLHS, outputs] = L_3D(consts, sim, laser, frame, mirror, outputs, toggles); % backward-propagating field

            % Update the boundary cells
            EplusNext(:,:,1) = EplusAtLHS;
            EminusNext(:,:,Nz) = EminusAtRHS;

            % Update the forward and backward propagating fields
            laser.Eplus = EplusNext;
            laser.Eminus = EminusNext;

            % Update clock
            sim.t = sim.t + sim.dt;

            % The seed reaches the RHS halfway through this round trip.
            if step == Nz && toggles.outputs_switch
                E_rhs = outputs.Gau_RHS_mirror;
                rhs_reflectance = outputs.R1(end);
            end

            Pplus = power_factor * reshape( ...
                sum(sum(abs(laser.Eplus).^2, 1), 2), 1, []);

            Pminus = power_factor * reshape( ...
                sum(sum(abs(laser.Eminus).^2, 1), 2), 1, []);

            outputs.t(end+1) = sim.t; % update time history, one per step
            outputs.Pplus(end+1,:) = Pplus;
            outputs.Pminus(end+1,:) = Pminus;
            outputs.U_cavity(end+1) = sim.dt*sum(Pplus + Pminus); % energy contained within the cavity

            if mod(step, frameStride) == 0
                Write_Cavity_Frame(video, fig, layout, sliceSurfaces, laser, ix, iy, sliceIndices, intensityFactor, sim.t);
            end

        end

        % At the end of every round-trip:

        E_end = laser.Eplus(:,:,1); % field after the round-trip is over

        if toggles.outputs_switch
            power_factor = 0.5*consts.c*consts.eps0*sim.dx*sim.dy;
            rt_loss = 1 - rhs_reflectance*outputs.R2(end); % round-trip loss
            outputs.loss_frac(end+1) = rt_loss;
            outputs.P_circ(end+1) = power_factor*sum(abs(E_end).^2, 'all'); % circulating power
            outputs.Gau_RHS = E_rhs;
            outputs.Gau_LHS = E_end;
            outputs.P_circ_RHS = power_factor*sum(abs(E_rhs).^2, 'all');
            outputs.P_circ_LHS = power_factor*sum(abs(E_end).^2, 'all');

            % Mode convergence check
            if check_convergence && isfinite(previous_loss) && isfinite(rt_loss) && any(E_start(:)) && any(E_end(:)) % only runs once some loss data is available
                [converged, state] = Check_Mode_Convergence(E_start, E_end, previous_loss, rt_loss, state, toggles, sim);
                if converged
                    break % finish simulation early
                end
            end

            previous_loss = rt_loss;
        end

    end

end