function Save_Variables(save_folder, varargin)

    % Ensure folder exists
    if ~exist(save_folder, 'dir')
        mkdir(save_folder);
    end

    % Loop over all input variables after the first
    for i = 1:numel(varargin)
        varname = inputname(i + 1); % get the name of the variable
        if isempty(varname)
            varname = sprintf('var%d', i);
        end

        filename = fullfile(save_folder, [varname '.mat']);
        value = varargin{i};
        save(filename, 'value');
        fprintf('Saved %s to %s\n', varname, filename);
    end
end
