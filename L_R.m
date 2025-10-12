% Laser travels from LHS mirror to RHS mirror

function [step, Z_traveled, Z_position, E, Es] = L_R(step, Z_traveled, Z_position, E, Es, save_interval, num_steps, dz, L, H, R, amask, tmask, grids)    

    for n = 1:num_steps
        
        step = step + 1;
        Z_traveled(step) = step*dz; % record total distance propagated so far
        
        % Record current position within the cavity, AFTER the step has
        % been taken
        if step == 1 % this is the first step
            Z_position(step) = -L/2 + dz;
        else
            Z_position(step) = Z_position(step-1) + dz;
        end

        % ---- manual pad (unchanged) ----
        padX = (grids.Nxp - grids.Nx)/2;
        padY = (grids.Nyp - grids.Ny)/2;
        Epad = zeros(grids.Nyp, grids.Nxp);
        Epad(padY+1 : padY+grids.Ny, padX+1 : padX+grids.Nx) = E;
        
        Epad = Epad .* grids.amask_pad;

        % ---- forward FFT (consistent shifts) ----
        % Note: compute FFT, then shift to center to match how grids.H was built
        Fe = fft2(Epad);          % raw FFT
        Fe = fftshift(Fe);        % center zero-frequency (now matches grids.H)

        % diagnostic 1: visualize magnitude in frequency space BEFORE propagation
        imagesc(log10(abs(Fe)+1e-12)); axis image; colorbar;
        title('log10 |F(kx,ky)| BEFORE H'); drawnow;
        
        % after multiplication
        Fe = Fe .* grids.H;   % or use Fe to overwrite
        imagesc(log10(abs(Fe)+1e-12)); axis image; colorbar;
        title('log10 |F(kx,ky)| AFTER H'); drawnow;
        
        % also print total power
        P_spatial = sum(abs(Epad(:)).^2);
        P_freq = sum(abs(Fe(:)).^2);
        fprintf('P spatial = %.6g, P freq = %.6g, ratio = %.3g\n', P_spatial, P_freq, P_spatial/P_freq);
        
        % ---- inverse FFT (undo shift correctly) ----
        Fe = ifftshift(Fe);       % move zero-frequency back to corner for ifft2
        Epad_out = ifft2(Fe);     % inverse transform (no extra fftshift)
        
        % ---- crop back to original and apply real-space masks ----
        Eout = Epad_out(padY+1 : padY+grids.Ny, padX+1 : padX+grids.Nx);
        E = Eout .* tmask .* amask;

        % % Propagation    
        % FE = fft2(E); % transform beam to frequency domain
        % FE = FE.*fftshift(H).*R(Z_position(step), Z_position(step)-dz); % propagate beam in frequency domain
        % E = ifft2(FE); % transform back to space domain
        % E = E.*tmask.*amask; % real operators, absorb energy at boundaries
        
        % Save E field snapshots
        % if mod(step, save_interval) == 0
        %     step_label = sprintf('step_%d', step);
        %     Es.(step_label) = E; % save intermediate field
        % end
    
    end

end