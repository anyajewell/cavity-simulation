function [consts, sim, laser, frame, mirror, outputs, gain_medium] = Load_Workspace(path, file)

    fullpath = fullfile(path, file);

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

    fprintf('Workspace successfully loaded from:\n%s\n', fullpath);

end