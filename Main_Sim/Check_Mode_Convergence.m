function [converged, state] = Check_Mode_Convergence(E_prev, E_curr, loss_prev, loss_curr, state)
    
    eps_rho = 1e-4; % overlap tolerance
    eps_lambda = 1e-4; % eigenvalue tolerance
    eps_loss = 1e-5; % loss tolerance
    N_consec = 3; % required consecutive passes

    if isempty(state)
        state.lambda_prev = [];
        state.counter = 0;
    end
    
    % Flatten fields for inner products
    E_prev_vec = E_prev(:);
    E_curr_vec = E_curr(:);
    
    % Test for field overlap
    num = abs(E_curr_vec' * E_prev_vec);
    den = norm(E_curr_vec) * norm(E_prev_vec);
    rho = num / den;
    
    overlap_ok = (1 - rho) < params.eps_rho;

    % Test for eigenvalue stability
    lambda = (E_curr_vec' * E_prev_vec) / (E_prev_vec' * E_prev_vec);
    
    if isempty(state.lambda_prev)
        lambda_ok = false; % can't evaluate yet
    else
        lambda_ok = abs(lambda - state.lambda_prev) / abs(state.lambda_prev) < params.eps_lambda;
    end

    state.lambda_prev = lambda;
    
    % Test for loss stability
    loss_ok = abs(loss_curr - loss_prev) < params.eps_loss;

    % Consecutive hit logic
    if overlap_ok && lambda_ok && loss_ok % all conditions satisfied
        state.counter = state.counter + 1;
    else
        state.counter = 0;
    end
    
    converged = state.counter >= N_consec;

end

