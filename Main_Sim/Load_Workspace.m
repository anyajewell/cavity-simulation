function [consts, sim, laser, frame, mirror, outputs, gain_medium, toggles] = Load_Workspace(fullpath)

    data = load(fullpath);

    % Validate required fields
    requiredFields = {'consts','sim','laser','frame','mirror','outputs','gain_medium','toggles'};
    
    for k = 1:length(requiredFields)
        if ~isfield(data, requiredFields{k})
            error('Missing required variable "%s" in workspace file.', requiredFields{k});
        end
    end

    % Assign outputs
    consts = data.consts;
    sim = data.sim;
    laser = data.laser;
    frame = data.frame;
    mirror = data.mirror;
    outputs = data.outputs;
    gain_medium = data.gain_medium;
    toggles = data.toggles;

    % Upgrade workspaces created before propagation used elapsed time as its step.
    if ~isfield(sim, 'direction')
        sim.direction = sign(sim.dz);
    end
    if ~isfield(sim, 'dt')
        sim.dt = abs(sim.dz)/consts.c;
    end
    if ~isfield(laser, 'time')
        laser.time = sim.t0;
    end
    if ~isfield(outputs, 'ts')
        outputs.ts = [];
    end

    fprintf('Workspace successfully loaded from:\n%s\n', fullpath);

end