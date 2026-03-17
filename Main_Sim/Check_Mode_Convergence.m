function [converged, state] = Check_Mode_Convergence(E_prev, E_curr, loss_prev, loss_curr, state, toggles, sim)
    
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
    
    overlap_ok = (1 - rho) < eps_rho;

    % Test for eigenvalue stability
    lambda = (E_curr_vec' * E_prev_vec) / (E_prev_vec' * E_prev_vec);
    
    if isempty(state.lambda_prev)
        lambda_ok = false; % can't evaluate yet
    else
        lambda_ok = abs(lambda - state.lambda_prev) / abs(state.lambda_prev) < eps_lambda;
    end

    state.lambda_prev = lambda;
    
    % Test for loss stability
    loss_ok = abs(loss_curr - loss_prev) < eps_loss;

    % Single-lobed / centered / radially-clean shape test
    mode_ok = Check_Single_Lobed_Mode(E_curr, sim.X, sim.Y);

    if strcmp(toggles.convergence_def, 'TEM00') % has to be lowest order mode
        if overlap_ok && lambda_ok && loss_ok && mode_ok
            state.counter = state.counter + 1;
        else
            state.counter = 0;
        end
    else
        if overlap_ok && lambda_ok && loss_ok
            state.counter = state.counter + 1;
        else
            state.counter = 0;
        end
    end
    
    converged = state.counter >= N_consec;

end

%% HELPER FUNCTION

function is_single_lobed = Check_Single_Lobed_Mode(E, X, Y)

    I = abs(E).^2;
    I = I / max(I(:));   % normalize intensity

    % Radius grid
    R = sqrt(X.^2 + Y.^2);

    % Characteristic grid spacing
    dx = mean(diff(X(1,:)));
    dy = mean(diff(Y(:,1)));
    dr = 0.5 * (abs(dx) + abs(dy));

    % Condition 1: Center brighter than surrounding annulus
    r_core  = 3 * dr;
    r_ring1 = 5 * dr;
    r_ring2 = 9 * dr;

    core_mask = R <= r_core;
    ring_mask = (R >= r_ring1) & (R <= r_ring2);

    I_core = mean(I(core_mask));
    I_ring = mean(I(ring_mask));

    center_dominant = I_core > I_ring;

    % Condition 2: No strong secondary peaks
    thresh = 0.1; % only look at peaks above 10% of max
    BW = imregionalmax(I) & (I > thresh);

    peak_vals = I(BW);
    peak_vals = sort(peak_vals, 'descend');

    if isempty(peak_vals)
        secondary_ok = false;
    elseif numel(peak_vals) == 1
        secondary_ok = true;
    else
        second_peak_ratio = peak_vals(2) / peak_vals(1);
        secondary_ok = second_peak_ratio < 0.3;
    end

    % Condition 3: Radial profile should be mostly monotone decreasing
    % This catches residual higher-order ring/halo structure.
    r_max = max(R(:));
    r_edges = 0:dr:r_max;
    r_centers = 0.5 * (r_edges(1:end-1) + r_edges(2:end));

    I_radial = zeros(size(r_centers));
    valid = false(size(r_centers));

    for k = 1:numel(r_centers)
        mask = (R >= r_edges(k)) & (R < r_edges(k+1));
        if any(mask(:))
            I_radial(k) = mean(I(mask));
            valid(k) = true;
        end
    end

    I_radial = I_radial(valid);

    % only assess where intensity is still meaningful
    keep = I_radial > 1e-3;
    I_radial = I_radial(keep);

    if numel(I_radial) < 3
        radial_ok = false;
    else
        dI = diff(I_radial);

        % measure how much the radial profile rises instead of falls
        upward_content = sum(max(dI, 0));
        total_variation = sum(abs(dI)) + eps;

        ring_metric = upward_content / total_variation;

        radial_ok = ring_metric < 0.05;  % profile can't rise much more than it falls
    end

    % Final heuristic
    is_single_lobed = center_dominant && secondary_ok && radial_ok;

end