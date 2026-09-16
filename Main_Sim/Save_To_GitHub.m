function gitOutDir = Save_To_GitHub(dir)
% Copy a local results folder into a local Git repository.

    % Make sure source folder exists
    if ~isfolder(dir)
        error('Source folder does not exist:\n%s', dir);
    end

    gitResultsDir = "C:\Users\jewellan\Documents\GitHub\cavity-simulation\Results\Dated_Data"; % general folder for dated data

    % Create Git Results directory if needed
    if ~isfolder(gitResultsDir)
        mkdir(gitResultsDir);
    end

    % Save
    [~, folderName] = fileparts(dir); % preserve the dated folder name
    gitOutDir = fullfile(gitResultsDir, folderName);
    [status, msg] = copyfile(dir, gitOutDir); % copy entire folder

    % Feedback Messages
    if ~status
        error('Failed to copy results to Git folder:\n%s', msg);
    end
    fprintf('Results copied to:\n%s\n', gitOutDir);

end