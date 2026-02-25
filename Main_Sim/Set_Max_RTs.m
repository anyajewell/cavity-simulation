function RTs = Set_Max_RTs(sim, consts, sampling_time)

    time_per_RT = sim.L/(consts.c); % amount of time it takes to travel one RT, [s]
    RTs = sampling_time / time_per_RT; % max number of RTs to propagate

end