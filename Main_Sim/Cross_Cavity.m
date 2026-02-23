% Start at either end and cross the cavity once. Output centers to compare
% with analytic ray optics solutions. This function is essentially half of
% Propagate_n_RTs.

function [laser, outputs, gain_medium] = Cross_Cavity(consts, sim, laser, frame, mirror, outputs, toggles, gain_medium)

    Gau_a = laser.Gau; % beam profile at the start

    if sim.dz > 0 % laser traveling L --> R
            [laser, outputs, sim] = R(consts, sim, laser, frame, mirror, outputs, toggles, gain_medium); % traveling right first
            outputs = Write_Video_Frame(sim, laser, toggles, outputs);

    else % laser traveling R --> L
        [laser, outputs, sim] = L(consts, sim, laser, frame, mirror, outputs, toggles, gain_medium); % traveling left first
        outputs = Write_Video_Frame(sim, laser, toggles, outputs);

    end
    
    if toggles.track_centers == true % fix counting
        outputs.centerx(end) = [];
        outputs.centery(end) = [];
    end

end