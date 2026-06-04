% Interface function between PCAC and laser propagation code

function [latest_P_circ, laser, outputs, gain_medium, sim] = Laser(dtheta_x, dtheta_y, consts, sim, laser, frame, mirror, outputs, toggles, gain_medium)

    % Update mirror misalignment from pointing error information
    mirror(1).dtheta_x = dtheta_x;
    mirror(1).dtheta_y = dtheta_y;

    [laser, outputs, gain_medium] = Propagate_n_RTs(consts, sim, laser, frame, mirror, outputs, toggles, gain_medium); % propagate until convergence
    
        if isfield(outputs, 'loss_frac') && ~isempty(outputs.loss_frac)
        
            P_circ_array = outputs.P_circ;
        
            % Get index of last nonzero entry
            idx = find(P_circ_array, 1, 'last'); % zeros are ignored automatically
        
            if ~isempty(idx)
                latest_P_circ = P_circ_array(idx); % deliver the converged loss fraction
            else
                % If array exists but is entirely zeros
                latest_P_circ = P_circ_array(end);
            end
        
        end

end