function Round_Trip(step, Z_traveled, Z_position, E, Es, save_interval, num_steps, dz, L, H)

    % Example of one round-trip
    
    % --------------------------------(1)----------------------------------- %
    
    % Beam enters cavity at +L/2 (RHS)
    
    % Propagation from +L/2 (RHS) to -L/2 (LHS)
    [step, Z_traveled, Z_position, E, Es] = R_L(step, Z_traveled, Z_position, E, Es, save_interval, num_steps, dz, L, H);
    
    % Beam interacts with mirror 2 (LHS) mirror at -L/2
    E = E.*cmask2.*rmask2;
    
    plot_intensity(E, consts, Z_traveled, Z_position, step, X, Y, x_circ2, y_circ2); % take a look at the beam
    
    % --------------------------------(2)----------------------------------- %
    
    % Propagation from -L/2 back to +L/2
    [step, Z_traveled, Z_position, E, Es] = L_R(step, Z_traveled, Z_position, E, Es, save_interval, num_steps, dz, L, H);
    
    % Beam interacts with mirror 1 (RHS) mirror at +L/2
    E = E.*cmask1.*rmask1;
    
    plot_intensity(E, consts, Z_traveled, Z_position, step, X, Y, x_circ2, y_circ2); % take a look at the beam

end