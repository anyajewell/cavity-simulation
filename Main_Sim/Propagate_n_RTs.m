function [laser, outputs, gain_medium] = Propagate_n_RTs(consts, sim, laser, frame, mirror, outputs, toggles, gain_medium)

    for a = 1:sim.RTs
        Gau_a = laser.Gau; % beam profile at the start of this round trip (RT)

        if sim.dz > 0 % laser traveling L --> R
            [laser, outputs, sim] = R(consts, sim, laser, frame, mirror, outputs, toggles, gain_medium); % traveling right first
            if toggles.outputs_switch == true && strcmp(toggles.videoplot_frequency, 'every mirror')
                outputs = Write_Video_Frame(sim, laser, toggles, outputs);
            end
            [laser, outputs, sim] = L(consts, sim, laser, frame, mirror, outputs, toggles, gain_medium); % traveling left
            if toggles.outputs_switch == true && strcmp(toggles.videoplot_frequency, 'every mirror')
                outputs = Write_Video_Frame(sim, laser, toggles, outputs);
            end
            
        else % laser traveling R --> L
            [laser, outputs, sim] = L(consts, sim, laser, frame, mirror, outputs, toggles, gain_medium); % traveling left first
            if toggles.outputs_switch == true && strcmp(toggles.videoplot_frequency, 'every mirror')
                outputs = Write_Video_Frame(sim, laser, toggles, outputs);
            end
            [laser, outputs, sim] = R(consts, sim, laser, frame, mirror, outputs, toggles, gain_medium); % traveling right
            if toggles.outputs_switch == true && strcmp(toggles.videoplot_frequency, 'every mirror')
                outputs = Write_Video_Frame(sim, laser, toggles, outputs);
            end

        end % One RT complete
        
        % Calculate and store RT loss
        loss_a =  1 - sum(abs(laser.Gau).^2,'all') / sum(abs(Gau_a).^2,'all'); % loss fraction this round-trip
        outputs.loss_frac(a) = loss_a;

    end

    if toggles.track_centers == true && toggles.outputs_switch == true % fix counting
        outputs.centerx(end) = [];
        outputs.centery(end) = [];
    end

end