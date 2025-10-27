% Laser travels from LHS mirror to RHS mirror

function [step, Z_traveled, Z_position, E, Es] = L_R(step, Z_traveled, Z_position, E, Es, save_interval, num_steps, dz, L, H, R, tmask, consts, N, Rdx_pixels, Omega, X, Y)    

    for n = 1:num_steps
        
        step = step + 1;
        Z_traveled(step) = step*dz; % record total distance propagated so far
        
        % Record current position within the cavity, AFTER the step has
        % been taken
        if step == 1 % this is the first step
            z0 = -L/2;
            z1 = z0 + dz;
            Z_position(step) = -L/2 + dz;
        else
            z0 = Z_position(step-1);
            z1 = z0 + dz;
            Z_position(step) = Z_position(step-1) + dz;
        end

        % Propagation    
        FE = fft2(E); % transform beam to frequency domain
        FE = FE.*fftshift(H); % propagate beam in frequency domain
        E = ifft2(FE); % transform back to space domain
        E = E.*tmask; % absorb energy at boundaries and apply tilting mask
        %alpha = (Omega/consts.c)*Z_position(step)*dz + (Omega/(2*consts.c))*dz^2;  % scalar shift in x
        alpha = (Omega/(2*consts.c))*(z1^2 - z0^2);  % scalar shift in x
        E = interp2(X, Y, E, X + alpha, Y, 'linear', 0); % shift the beam, rotational shearing
        %E = interp2(E, (1:N) - Rdx_pixels(Z_position(step)), (1:N)', 'linear', 0); % shift the beam, rotational shearing
        
        % % Save E field snapshots
        % if mod(step, save_interval) == 0
        %     step_label = sprintf('step_%d', step);
        %     Es.(step_label) = E; % save intermediate field
        % end
        % 
    end

end