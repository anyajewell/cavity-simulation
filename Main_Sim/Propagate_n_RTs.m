function [laser, outputs, gain_medium] = Propagate_n_RTs(consts, sim, laser, frame, mirror, outputs, toggles, gain_medium)

    for a = 1:sim.RTs
        Gau_a = laser.Gau; % beam profile at the start of this round trip

        if sim.dz > 0 % laser traveling L --> R
            [laser, outputs, sim] = R(consts, sim, laser, frame, mirror, outputs, toggles, gain_medium, a); % traveling right first
            outputs = Write_Video_Frame(sim, laser, toggles, outputs);
            [laser, outputs, sim] = L(consts, sim, laser, frame, mirror, outputs, toggles, gain_medium, a); % traveling left
            outputs = Write_Video_Frame(sim, laser, toggles, outputs);
            
        else % laser traveling R --> L
            [laser, outputs, sim] = L(consts, sim, laser, frame, mirror, outputs, toggles, gain_medium, a); % traveling left first
            outputs = Write_Video_Frame(sim, laser, toggles, outputs);
            [laser, outputs, sim] = R(consts, sim, laser, frame, mirror, outputs, toggles, gain_medium, a); % traveling right
            outputs = Write_Video_Frame(sim, laser, toggles, outputs);

        end
        
        % Calculate and store RT loss
        loss_a =  1 - sum(abs(laser.Gau).^2,'all') / sum(abs(Gau_a).^2,'all'); % loss fraction this round-trip
        outputs.loss_frac(a) = loss_a;
    end

    if toggles.track_centers == true % fix counting
        outputs.centerx(end) = [];
        outputs.centery(end) = [];
    end

end