function [Gau, loss_frac] = Propagate_n_RTs(RTs, Gau, Nz, Omega, accel, dt, c, Ld, dx, dz, x, y, z, X, Y, k0, v, centerx, centery, cmask1, cmask2, rmask1, rmask2, L, Zmax)

    for a = 1:RTs
        Gau_a = Gau; % beam profile at the start of this round trip

        % L --> R
        for i = 1:Nz
            pos = z(i); % current position
            Omega = Omega; % opportunity to read/update Omega from other code
            accel = accel; % opportunity to read/update accel from other code
            if pos + dz > L/2 
                dz1 = Zmax-pos; 
                if dz1 == 0 % beam already at mirror
                    break
                else % propagate by fractional step to L (?)
                [Gau, centerx, centery, v] = Prop(Gau, Omega, accel, dt, c, Ld, dx, dz1, x, y, z, X, Y, k0, v, centerx, centery, i); % propagate fractional step to mirror
                end
            else % propagate by dz, as normal
                [Gau, centerx, centery, v] = Prop(Gau, Omega, accel, dt, c, Ld, dx, dz, x, y, z, X, Y, k0, v, centerx, centery, i); % propagation loop
            end
        end
    
        % Interact with mirror (RHS)
        Gau = Gau.*cmask1.*rmask1; % clip and shape the beam
        [Zmax, dz, z] = Turn_Around(Zmax, dz, Nz);
    
        % R --> L
        for i = 1:Nz
            pos = z(i); % current position
            Omega = Omega; % opportunity to read/update Omega from other code
            accel = accel; % opportunity to read/update accel from other code
            if pos + dz < -L/2 
                dz1 = Zmax-pos; 
                if dz1 == 0 % beam already at mirror
                    break
                else % propagate by fractional step to L (?)
                [Gau, centerx, centery, v] = Prop(Gau, Omega, accel, dt, c, Ld, dx, dz1, x, y, z, X, Y, k0, v, centerx, centery, i); % propagate fractional step to mirror
                end
            else % propagate by dz, as normal
                [Gau, centerx, centery, v] = Prop(Gau, Omega, accel, dt, c, Ld, dx, dz, x, y, z, X, Y, k0, v, centerx, centery, i); % propagation loop
            end
        end
    
        % Interact with mirror (LHS)
        Gau = Gau.*cmask2.*rmask2; % clip and shape the beam
        [Zmax, dz, z] = Turn_Around(Zmax, dz, Nz);
    
        % Calculate and store loss
        loss_a =  1 - sum(abs(Gau),'all') / sum(abs(Gau_a).^2,'all'); % loss this time
        loss_frac(a) = loss_a;
    end

end