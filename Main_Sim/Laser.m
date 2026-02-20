% Interface function between PCAC and laser propagation code

function [loss_frac, laser, outputs, gain_medium] = Laser(dtheta_x, dtheta_y, sampling_time, consts, sim, laser, frame, mirror, outputs, toggles, gain_medium)

    mirror(1).dtheta_x = dtheta_x;
    mirror(1).dtheta_y = dtheta_y;
    time_per_RT = sim.L/(consts.c); % amount of time it takes to travel one RT, [s]
    sim.RTs = sampling_time / time_per_RT; % number of RTs to propagate
    [laser, outputs, gain_medium] = Propagate_n_RTs(consts, sim, laser, frame, mirror, outputs, toggles, gain_medium);
    loss_frac = outputs.loss_frac(end);

end