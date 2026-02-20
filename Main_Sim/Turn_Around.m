function [sim] = Turn_Around(sim)

    if pos == mirror(1).loc % laser at RHS of cavity
        sim.dz = -sim.dz;
        sim.z = linspace(mirror(1).loc,mirror(2).loc,sim.Nz);
    else % laser at LHS of cavity
        sim.dz = -sim.dz;
        sim.z = linspace(mirror(2).loc,mirror(1).loc,sim.Nz);
    end

end