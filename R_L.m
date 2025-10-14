% Laser travels from RHS mirror to LHS mirror

function [step, Z_traveled, Z_position, E, Es] = R_L(step, Z_traveled, Z_position, E, Es, save_interval, num_steps, dz, L, H, R, amask, consts)

    for n = 1:num_steps
        
        step = step + 1;
        Z_traveled(step) = step*dz; % record total distance propagated so far

        % Record current position within the cavity, AFTER the step has
        % been taken
        if step == 1 % this is the first step
            Z_position(step) = L/2 - dz;
        else
            Z_position(step) = Z_position(step-1) - dz;
        end

        % Propagation
        FE = fft2(E); % transform beam to frequency domain
        FE = FE.*fftshift(H).*R(Z_position(step), Z_position(step)-dz); % propagate beam in frequency domain 
        E = ifft2(FE); % transform back to space domain 
        E = E.*amask.*tmask; % absorb energy at boundaries and apply tilting mask

        % % Save E field snapshots
        % if mod(step, save_interval) == 0
        %     step_label = sprintf('step_%d', step);
        %     Es.(step_label) = E; % save intermediate field
        % end
    
    end
end