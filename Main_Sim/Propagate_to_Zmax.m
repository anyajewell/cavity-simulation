function [laser, outputs] = Propagate_to_Zmax(Z0, Zmax, consts, sim, laser, frame, outputs, toggles)
    
    dz = Zmax/sim.Nz;
    z = linspace(Z0, Zmax, sim.Nz);
    sim.dz = dz; sim.z = z;
    laser.pos = Z0;

    while laser.pos < Zmax
        [laser, outputs] = Prop(consts, sim, laser, frame, outputs, toggles, sim.dz); % propagation loop
        laser.pos = laser.pos + sim.dz; % update laser position   
        outputs = Write_Video_Frame(sim, laser, toggles, outputs);
    end

    outputs.w_num = sqrt(2*sum(sum((sim.X.^2 + sim.Y.^2).*abs(laser.Gau).^2))/sum(sum(abs(laser.Gau).^2)));
    outputs.P_num = sum(sum(abs(laser.Gau).^2))*sim.dx^2;
    R = sqrt(sim.X.^2 + sim.Y.^2);
    [~,idx] = min(abs(R(:)-outputs.w_num)); % picks out the location of beam waist
    I_ratio = abs(laser.Gau(idx)).^2 / max(abs(laser.Gau(:)).^2); % should be 0.1353

end
