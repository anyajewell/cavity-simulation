function latest_loss = Get_Latest_Loss_Fraction(matfile)

    if ~isfile(matfile)
        error('File not found: %s', matfile);
    end

    data = load(matfile);
    vars = fieldnames(data);

    latest_loss = [];

    for k = 1:length(vars)
        varname = vars{k};
        value = data.(varname);

        % Check if variable is a struct with loss_frac field
        if isstruct(value) && isfield(value, 'loss_frac')
            
            loss_array = value.loss_frac;

            if ~isempty(loss_array)
                % Remove trailing zeros if preallocated
                nonzero_idx = find(loss_array ~= 0);
                
                if ~isempty(nonzero_idx)
                    latest_loss = loss_array(nonzero_idx(end));
                else
                    latest_loss = loss_array(end);
                end
                
                return
            end
        end
    end

    error('No loss_frac field found in any struct within %s.', matfile);

end