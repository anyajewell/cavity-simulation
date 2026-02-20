function Save_Workspace(consts, sim, laser, frame, mirror, outputs, gain_medium)

    % Ensure folder exists
    if ~exist(outputs.saveFolder, 'dir')
        mkdir(outputs.saveFolder);
    end

    % Timestamp for uniqueness
    timestamp = datestr(now,'yyyy_mm_dd_HHMMSS');

    % File name
    filename = fullfile(outputs.saveFolder, ['workspace_' timestamp '.mat']);

    % Save all major structs
    save(filename, ...
        'consts', ...
        'sim', ...
        'laser', ...
        'frame', ...
        'mirror', ...
        'outputs', ...
        'gain_medium', ...
        '-v7.3');   % v7.3 handles large arrays safely

    fprintf('Workspace saved to:\n%s\n', filename);

end
