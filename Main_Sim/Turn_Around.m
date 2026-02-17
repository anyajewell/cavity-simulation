function [sim] = Turn_Around(sim)

    % Turn around
    sim.Zmax = -sim.Zmax;
    sim.Z0 = -sim.Zmax;
    sim.dz = -sim.dz;
    sim.z = linspace(sim.Z0,sim.Zmax,sim.Nz);

end