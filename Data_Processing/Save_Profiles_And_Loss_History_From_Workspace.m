clear;
clc;
close all;

% Base directory
base_folder = 'C:\Users\jewellan\Documents\GitHub\cavity-simulation\Results\SciTech_Data\Fresnel_Omega_Study';

% Sweep parameters
NF_list = 0.1:0.1:1.0;
Omega_list = 0.02:0.01:0.1;

total_cases = numel(NF_list)*numel(Omega_list);

for a = 1:length(Omega_list)

    Omega = Omega_list(a);

    % Omega folder
    root_folder = fullfile(base_folder, ...
        sprintf('Fresnel_Sweep_Omega=%.3f', Omega));

    for i = 1:length(NF_list)

        N_F = NF_list(i);
        case_num = (a-1)*numel(NF_list) + i;

        % Fresnel number folder
        folder_name = sprintf('NF_%0.1f', N_F);
        folder_name = strrep(folder_name, '.', 'p');

        case_folder = fullfile(root_folder, folder_name);
        workspace_file = fullfile(case_folder, 'workspace.mat');

        fprintf('\n====================================\n');
        fprintf('Processing case %d/%d\n', case_num, total_cases);
        fprintf('Omega = %.3f, NF = %.1f\n', Omega, N_F);
        fprintf('====================================\n');

        % Check workspace exists
        if ~isfile(workspace_file)
            warning('Workspace missing: %s', workspace_file);
            continue;
        end

        % Load only required variables
        data = load(workspace_file, 'sim', 'outputs');

        sim = data.sim;
        outputs = data.outputs;

        close all;

        %% RHS laser profile
        if isfield(outputs, 'Gau_RHS')

            Plot_Laser_Profile(sim, outputs.Gau_RHS);

            exportgraphics(gcf, ...
                fullfile(case_folder, 'Gau_RHS.png'), ...
                'Resolution', 300);

            close(gcf);
        end

        %% LHS laser profile
        if isfield(outputs, 'Gau_LHS')

            Plot_Laser_Profile(sim, outputs.Gau_LHS);

            exportgraphics(gcf, ...
                fullfile(case_folder, 'Gau_LHS.png'), ...
                'Resolution', 300);

            close(gcf);
        end

        %% Round-trip loss
        if isfield(outputs, 'loss_frac') && ~isempty(outputs.loss_frac)

            Plot_RT_Loss(outputs.loss_frac, length(outputs.loss_frac));

            exportgraphics(gcf, ...
                fullfile(case_folder, 'RT_Loss.png'), ...
                'Resolution', 300);

            close(gcf);
        end

        fprintf('Plots saved to:\n%s\n', case_folder);

        % Clear loaded simulation data
        clear sim outputs data

    end
end

fprintf('\nAll cases processed.\n');