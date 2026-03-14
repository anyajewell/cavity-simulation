function gauss_debug = Analytic_Gaussian(sim, laser, mirror)
    
    % Analytic Gaussian-beam comparison for: propagate LHS to RHS, apply
    % RHS mirror, propagate RHS to LHS, apply LHS mirror
   
    % q0: the initial beam parameter
    % E0: the initial on-axis amplitude
    % At the initial plane: E = E0 * exp(-i*k*r^2/(2*q0))
    % Free-space propagation: q_out = q_in + L
    % Mirror reflection: 1/q_out = 1/q_in - 2/Rm or q_out = q_in / (1 - 2*q_in/Rm)

    E0 = 1;
    q0 = 1i*pi*laser.w0^2/laser.Ld; % initial q-parameter
    k = 2*pi/laser.Ld;
    r2 = sim.X.^2 + sim.Y.^2;
    gauss_debug.lambda = laser.Ld;
    gauss_debug.L = sim.L;
    gauss_debug.R1 = mirror(1).Rc;
    gauss_debug.R2 = mirror(2).Rc;

    % Evolution of q-parameter
    q_LHS_0 = q0; % initial Gaussian at LHS
    q_RHS_before = q_LHS_0 + sim.L; % propagate LHS -> RHS
    q_RHS_after = q_RHS_before ./ (1 - 2*q_RHS_before/mirror(2).Rc); % apply RHS mirror
    q_LHS_return = q_RHS_after + sim.L; % propagate RHS -> LHS
    q_LHS_after_LHSmirror = q_LHS_return ./ (1 - 2*q_LHS_return/mirror(1).Rc); % apply LHS mirror to complete round trip

    % Store q values
    gauss_debug.q_LHS_0 = q_LHS_0;
    gauss_debug.q_RHS_before = q_RHS_before;
    gauss_debug.q_RHS_after = q_RHS_after;
    gauss_debug.q_LHS_return = q_LHS_return;
    gauss_debug.q_LHS_after_LHSmirror = q_LHS_after_LHSmirror;

    % Convert q -> w and R
    [gauss_debug.w_LHS_0, gauss_debug.R_LHS_0] = q_to_wR(q_LHS_0, laser.Ld);
    [gauss_debug.w_RHS_before, gauss_debug.R_RHS_before] = q_to_wR(q_RHS_before, laser.Ld);
    [gauss_debug.w_RHS_after, gauss_debug.R_RHS_after] = q_to_wR(q_RHS_after, laser.Ld);
    [gauss_debug.w_LHS_return, gauss_debug.R_LHS_return] = q_to_wR(q_LHS_return, laser.Ld);
    [gauss_debug.w_LHS_after_LHSmirror, gauss_debug.R_LHS_after_LHSmirror] = ...
        q_to_wR(q_LHS_after_LHSmirror, laser.Ld);

    % Analytic complex field profiles
    gauss_debug.E_LHS_0 = gaussian_field_from_q(r2, k, q_LHS_0, q0, E0);
    gauss_debug.E_RHS_before = gaussian_field_from_q(r2, k, q_RHS_before, q0, E0);
    gauss_debug.E_RHS_after = gaussian_field_from_q(r2, k, q_RHS_after, q0, E0);
    gauss_debug.E_LHS_return = gaussian_field_from_q(r2, k, q_LHS_return, q0, E0);
    gauss_debug.E_LHS_after_LHSmirror = gaussian_field_from_q(r2, k, q_LHS_after_LHSmirror, q0, E0);

    % Print compact summary
    fprintf('\n=== Analytic Gaussian Debug (with field) ===\n');
    fprintf('Initial LHS:\n');
    fprintf('  q = %.6e %+.6ei\n', real(q_LHS_0), imag(q_LHS_0));
    fprintf('  w = %.6e m, R = %.6e m\n', gauss_debug.w_LHS_0, gauss_debug.R_LHS_0);

    fprintf('\nAfter LHS -> RHS propagation:\n');
    fprintf('  q = %.6e %+.6ei\n', real(q_RHS_before), imag(q_RHS_before));
    fprintf('  w = %.6e m, R = %.6e m\n', gauss_debug.w_RHS_before, gauss_debug.R_RHS_before);

    fprintf('\nAfter RHS mirror:\n');
    fprintf('  q = %.6e %+.6ei\n', real(q_RHS_after), imag(q_RHS_after));
    fprintf('  w = %.6e m, R = %.6e m\n', gauss_debug.w_RHS_after, gauss_debug.R_RHS_after);

    fprintf('\nAfter RHS -> LHS propagation:\n');
    fprintf('  q = %.6e %+.6ei\n', real(q_LHS_return), imag(q_LHS_return));
    fprintf('  w = %.6e m, R = %.6e m\n', gauss_debug.w_LHS_return, gauss_debug.R_LHS_return);

    fprintf('\nAfter LHS mirror (full RT):\n');
    fprintf('  q = %.6e %+.6ei\n', real(q_LHS_after_LHSmirror), imag(q_LHS_after_LHSmirror));
    fprintf('  w = %.6e m, R = %.6e m\n', gauss_debug.w_LHS_after_LHSmirror, gauss_debug.R_LHS_after_LHSmirror);
    fprintf('============================================\n\n');
end

% HELPER FUNCTIONS

function E = gaussian_field_from_q(r2, k, q, q_ref, E0)
    E = E0 .* (q_ref./q) .* exp(-1i*k*r2./(2*q)); % returns the complex Gaussian field profile at a plane with beam parameter q
end

function [w, R] = q_to_wR(q, lambda)
    
    invq = 1./q;
    w = sqrt(-lambda ./ (pi * imag(invq))); % spot size from imaginary part of 1/q

    % Convert complex q parameter into spot size and ROC
    if abs(real(invq)) < 1e-14
        R = Inf;
    else
        R = 1 ./ real(invq);
    end

end