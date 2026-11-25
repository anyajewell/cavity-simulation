clear;
clc;
close all;

% Base directory
base_folder = 'C:\Users\jewellan\Documents\GitHub\cavity-simulation\Results\SciTech_Data\Fresnel_Omega_Study';

% Sweep parameters
NF_list = 0.1:0.1:1.0;
Omega_list = 0:0.001:0.01;
Omega_list2 = 0.02:0.01:0.1;

Omega_list = [Omega_list, Omega_list2];

% Preallocate loss matrix
% Rows = Omega values, Columns = Fresnel numbers
loss_matrix = NaN(length(Omega_list), length(NF_list));

%% Extract converged loss from each case

for a = 1:length(Omega_list)

    Omega = Omega_list(a);

    root_folder = fullfile(base_folder, ...
        sprintf('Fresnel_Sweep_Omega=%.3f', Omega));

    for i = 1:length(NF_list)

        N_F = NF_list(i);

        folder_name = sprintf('NF_%0.1f', N_F);
        folder_name = strrep(folder_name, '.', 'p');

        case_folder = fullfile(root_folder, folder_name);
        workspace_file = fullfile(case_folder, 'workspace.mat');

        if ~isfile(workspace_file)
            warning('Missing workspace: %s', workspace_file);
            continue;
        end

        % Load outputs only
        data = load(workspace_file, 'outputs');

        if isfield(data.outputs, 'loss_frac') && ...
                ~isempty(data.outputs.loss_frac)

            % Final round-trip loss, converted to percent
            loss_matrix(a,i) = data.outputs.loss_frac(end)*100;

        else
            warning('Missing loss data: %s', workspace_file);
        end

    end
end

%% Plot loss versus Fresnel number

figure('Color', 'w', 'Position', [100 100 900 600]);
hold on;

% Continuous colormap for Omega
cmap = turbo(256);
colormap(cmap);

Omega_min = min(Omega_list);
Omega_max = max(Omega_list);

% Map each Omega to its corresponding color
for a = 1:length(Omega_list)

    Omega = Omega_list(a);

    % Normalize Omega to [0, 1]
    normalized_Omega = (Omega - Omega_min) / (Omega_max - Omega_min);

    % Convert to colormap index
    color_idx = round(normalized_Omega * 255) + 1;
    color = cmap(color_idx, :);

    % Plot
    semilogy(NF_list, loss_matrix(a,:), ...
        '-o', ...
        'Color', color, ...
        'LineWidth', 1.8, ...
        'MarkerSize', 4, ...
        'MarkerFaceColor', color);

end

% Colorbar

clim([Omega_min Omega_max]);

cb = colorbar;
cb.Label.String = '\Omega [rad/s]';
cb.Label.FontSize = 12;
cb.FontSize = 11;

% Formatting

xlabel('N_F', 'FontSize', 12);
ylabel('Loss / round-trip (%)', 'FontSize', 12);

title({'Fresnel number versus loss per round-trip', ...
       'for varying rotation rates'}, ...
       'FontSize', 13);

grid on;
grid minor;
box on;

set(gca, 'FontSize', 11, 'YScale', 'log');

xlim([min(NF_list), max(NF_list)]);

hold off;

%% Save combined figure

exportgraphics(gcf, ...
    fullfile(base_folder, 'Fresnel_Sweep_Converged_Loss.png'), ...
    'Resolution', 300);

savefig(gcf, ...
    fullfile(base_folder, 'Fresnel_Sweep_Converged_Loss.fig'));

% Save numerical data for later analysis
save(fullfile(base_folder, 'Fresnel_Sweep_Loss_Data.mat'), ...
    'NF_list', 'Omega_list', 'loss_matrix');

fprintf('\nCombined loss plot and data saved successfully.\n');