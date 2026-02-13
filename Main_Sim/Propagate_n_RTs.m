function [Gau, loss_frac, centerx, centery] = Propagate_n_RTs(RTs, Gau, Nz, Omega, accel, dt, c, Ld, dx, dz, x, y, z, X, Y, k0, v, centerx, centery, cmask1, cmask2, rmask1, rmask2, L, Zmax, track_centers)

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
                    if track_centers == true
                        centerx(end+1) = trapz(trapz(X.*abs(Gau).^2))/trapz(trapz(abs(Gau).^2)); % track center x
                        centery(end+1) = trapz(trapz(Y.*abs(Gau).^2))/trapz(trapz(abs(Gau).^2)); % track center y
                    end
                    break
                else % propagate by fractional step to L
                [Gau, centerx, centery, v] = Prop(Gau, Omega, accel, dt, c, Ld, dx, dz1, x, y, z, X, Y, k0, v, centerx, centery, i, track_centers); % propagate fractional step to mirror
                end
            else % propagate by dz, as normal
                [Gau, centerx, centery, v] = Prop(Gau, Omega, accel, dt, c, Ld, dx, dz, x, y, z, X, Y, k0, v, centerx, centery, i, track_centers); % propagation loop
            end
        end
    
        % Interact with mirror (RHS)
        theta_x1 = 0; theta_y1 = 0; % query mirror misalignment
        phi_tilt1 = 2*k0*(theta_x1*X + theta_y1*Y);
        tmask1 = exp(1i*phi_tilt1); % tilt mask
        Gau = Gau.*cmask1.*rmask1.*tmask1; % clip and shape the beam
        [Zmax, dz, z] = Turn_Around(Zmax, dz, Nz);

        imagesc(x,y,abs(Gau)); axis([-0.5 0.5 -0.5 0.5]); axis square; xlabel('x [m]'); ylabel('y [m]'); hold on; 
        if track_centers == true
            plot(centerx(end),centery(end),'ro'); hold off;
        end
        
        frame = getframe(gcf); display(z(i));
        writeVideo(v,frame);
    
        % R --> L
        for i = 1:Nz
            pos = z(i); % current position
            Omega = Omega; % opportunity to read/update Omega from other code
            accel = accel; % opportunity to read/update accel from other code
            if pos + dz < -L/2 
                dz1 = Zmax-pos; 
                if dz1 == 0 % beam already at mirror
                    if track_centers == true
                        centerx(end+1) = trapz(trapz(X.*abs(Gau).^2))/trapz(trapz(abs(Gau).^2)); % track center x
                        centery(end+1) = trapz(trapz(Y.*abs(Gau).^2))/trapz(trapz(abs(Gau).^2)); % track center y
                    end
                    break
                else % propagate by fractional step to L
                [Gau, centerx, centery, v] = Prop(Gau, Omega, accel, dt, c, Ld, dx, dz1, x, y, z, X, Y, k0, v, centerx, centery, i, track_centers); % propagate fractional step to mirror
                end
            else % propagate by dz, as normal
                [Gau, centerx, centery, v] = Prop(Gau, Omega, accel, dt, c, Ld, dx, dz, x, y, z, X, Y, k0, v, centerx, centery, i, track_centers); % propagation loop
            end
        end
    
        % Interact with mirror (LHS)
        theta_x2 = 0; theta_y2 = 0; % query mirror misalignment
        if a == 30
            theta_x2 = 0.5e-6; theta_y2 = 0.5e-6;
        end
        phi_tilt2 = 2*k0*(theta_x2*X + theta_y2*Y);
        tmask2 = exp(1i*phi_tilt2); % tilt mask
        Gau = Gau.*cmask2.*rmask2.*tmask2; % clip and shape the beam
        [Zmax, dz, z] = Turn_Around(Zmax, dz, Nz);

        % Visualization
        imagesc(x,y,abs(Gau)); axis([-0.5 0.5 -0.5 0.5]); axis square; xlabel('x [m]'); ylabel('y [m]');  
        set(gcf, 'Color', 'w'); set(gca, 'Color', 'w'); hold on;
        if track_centers == true
            plot(centerx(end),centery(end),'ro'); hold off;
        end
        frame = getframe(gcf); display(z(i));
        writeVideo(v,frame);
    
        % Calculate and store loss
        loss_a =  1 - sum(abs(Gau).^2,'all') / sum(abs(Gau_a).^2,'all'); % loss this time
        loss_frac(a) = loss_a;
    end

    if track_centers == true % fix counting
        centerx(end) = [];
        centery(end) = [];
    end

end