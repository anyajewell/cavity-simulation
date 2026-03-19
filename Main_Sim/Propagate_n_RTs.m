function [laser, outputs, gain_medium] = Propagate_n_RTs(consts, sim, laser, frame, mirror, outputs, toggles, gain_medium)

    state = []; % initialize for convergence checks

    for a = 1:sim.RTs

        laser.Gau_a = laser.Gau; % beam profile at the start of this round trip (RT)

        if sim.dz > 0 % laser traveling L --> R
            [laser, outputs, sim] = R(consts, sim, laser, frame, mirror, outputs, toggles, gain_medium); % traveling right first
            Gau_RHS = laser.Gau; % profile at mirror 1
            if a == 1 && toggles.resize_grid == true
                [sim, laser, mirror, gain_medium] = Decrease_To_Grid0_Centroid(sim, laser, mirror, outputs, toggles);
            end
            if toggles.outputs_switch == true && strcmp(toggles.videoplot_frequency, 'every mirror')
                outputs = Write_Video_Frame(sim, laser, toggles, outputs);
            end
            [laser, outputs, sim] = L(consts, sim, laser, frame, mirror, outputs, toggles, gain_medium); % traveling left
            Gau_LHS = laser.Gau; % profile at mirror 2
            if toggles.outputs_switch == true && strcmp(toggles.videoplot_frequency, 'every mirror')
                outputs = Write_Video_Frame(sim, laser, toggles, outputs);
            end
            
        else % laser traveling R --> L
            [laser, outputs, sim] = L(consts, sim, laser, frame, mirror, outputs, toggles, gain_medium); % traveling left first
            Gau_LHS = laser.Gau; % profile at mirror 2
            if a == 1 && toggles.resize_grid == true
                [sim, laser, mirror, gain_medium] = Decrease_To_Grid0_Centroid(sim, laser, mirror, outputs, toggles);
            end
            if toggles.outputs_switch == true && strcmp(toggles.videoplot_frequency, 'every mirror')
                outputs = Write_Video_Frame(sim, laser, toggles, outputs);
            end
            [laser, outputs, sim] = R(consts, sim, laser, frame, mirror, outputs, toggles, gain_medium); % traveling right
            Gau_RHS = laser.Gau; % profile at mirror 1
            if toggles.outputs_switch == true && strcmp(toggles.videoplot_frequency, 'every mirror')
                outputs = Write_Video_Frame(sim, laser, toggles, outputs);
            end

        end % One RT complete
        
        % Calculate and store RT loss
        if isfield(laser,'P_ref') && numel(laser.P_ref) >= 2 
            loss_a = 1 - (1-outputs.loss2(end))*(1-outputs.loss1(end)); % round-trip loss
            outputs.loss_frac(a) = loss_a;
        end

        % Mode convergence check
        if a >= 2 && strcmp(toggles.finish_line, 'convergence') % only run once some loss data is available
            [converged, state] = Check_Mode_Convergence(laser.Gau_a, laser.Gau, outputs.loss_frac(a-1), outputs.loss_frac(a), state, toggles, sim);
            if converged == true || a == sim.RTs % mode convergence detected or simulation over
                outputs.Gau_RHS = Gau_RHS;
                outputs.Gau_LHS = Gau_LHS;
                break % finish propagation early (time skip)
            end
        end

    % if mod(a,20)==0
    %     filename = fullfile(outputs.saveFolder, 'checkpoint.mat');
    %     save(filename, 'outputs', 'laser', 'sim');
    % end

    end

    if toggles.track_centers == true && toggles.outputs_switch == true % fix counting
        outputs.centerx(end) = [];
        outputs.centery(end) = [];
    end

end