function latest_P_circ = Get_Latest_P_Circ(matfile)

    if ~isfile(matfile)
        error('File not found: %s', matfile);
    end

    data = load(matfile);
    vars = fieldnames(data);

    latest_P_circ = [];

    for k = 1:length(vars)
        varname = vars{k};
        value = data.(varname);

        % Check if variable is a struct with P_circ field
        if isstruct(value) && isfield(value, 'P_circ')
            
            P_circ_array = value.loss_frac;

            if ~isempty(P_circ_array)
                % Remove trailing zeros if preallocated
                nonzero_idx = find(P_circ_array ~= 0);
                
                if ~isempty(nonzero_idx)
                    latest_P_circ = P_circ_array(nonzero_idx(end));
                else
                    latest_P_circ = P_circ_array(end);
                end
                
                return
            end
        end
    end

    error('No P_circ field found in any struct within %s.', matfile);

end