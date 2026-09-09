function loss_array = Get_Final_Losses_Array(parent_folder)

    % Find all NF folders
    folders = dir(fullfile(parent_folder, 'NF_*'));
    folders = folders([folders.isdir]);

    loss_array = zeros(1, length(folders));

    for i = 1:length(folders)

        folder_path = fullfile(parent_folder, folders(i).name);

        % Find MAT files in this folder
        files = dir(fullfile(folder_path, '*.mat'));

        % Use the most recently modified MAT file
        [~, idx] = max([files.datenum]);
        file_path = fullfile(folder_path, files(idx).name);

        % Load file
        data = load(file_path);

        % Grab final loss value
        loss_array(i) = data.outputs.loss_frac(end);

    end

end