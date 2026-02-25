% Interface function between PCAC and laser propagation code

function [latest_loss, laser, outputs, gain_medium, sim] = Laser(dtheta_x, dtheta_y, consts, sim, laser, frame, mirror, outputs, toggles, gain_medium)

    % Update mirror misalignment from pointing error information
    mirror(1).dtheta_x = dtheta_x;
    mirror(1).dtheta_y = dtheta_y;

    outputs.loss_frac = zeros(1, sim.RTs); % re-allocate loss-frac array to be filled

    [laser, outputs, gain_medium] = Propagate_n_RTs(consts, sim, laser, frame, mirror, outputs, toggles, gain_medium); % propagate until convergence

    % Deliver the converged loss fraction
    latest_loss = [];
    
        if isfield(outputs, 'loss_frac') && ~isempty(outputs.loss_frac)
        
            loss_array = outputs.loss_frac;
        
            % Get index of last nonzero entry
            idx = find(loss_array, 1, 'last'); % zeros are ignored automatically
        
            if ~isempty(idx)
                latest_loss = loss_array(idx);
            else
                % If array exists but is entirely zeros
                latest_loss = loss_array(end);
            end
        
        end

end