% Laser travels from LHS mirror to RHS mirror

function [step, Z_traveled, Z_position, E, Es] = L_R(step, Z_traveled, Z_position, E, Es, save_interval, num_steps, dz, L, H)    

    for n = 1:num_steps
        
        step = step + 1;
        Z_traveled(step) = step*dz; % record total distance propagated so far
        
        % record current position within the cavity, AFTER the step has
        % been taken
        if step == 1 % this is the first step
            Z_position(step) = -L/2 + dz;
        else
            Z_position(step) = Z_position(step-1) + dz;
        end

        % propagation    
        FE = fft2(E); % transform beam to frequency domain
        FE = FE.*fftshift(H); % propagate beam in frequency domain 
        E = ifft2(FE); % transform back to space domain 
        
        % save E field snapshots
        if mod(step, save_interval) == 0
            step_label = sprintf('step_%d', step);
            Es.(step_label) = E; % save intermediate field
        end
    
    end

end