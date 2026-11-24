clear;
clc;
close all;

%% Directories and sweep parameters

base_folder = ...
    'C:\Users\jewellan\Documents\GitHub\cavity-simulation\Results\SciTech_Data\Fresnel_Omega_Study';

save_folder = fullfile(base_folder, 'Mode_Comparison_Figures');

if ~exist(save_folder, 'dir')
    mkdir(save_folder);
end

NF_list = 0.1:0.1:1.0;
Omega_list = 0.001:0.001:0.01;

nNF = numel(NF_list);
nOmega = numel(Omega_list);

%% Load all mode profiles

Gau_RHS = cell(nOmega, nNF);
Gau_LHS = cell(nOmega, nNF);

sim_ref = [];

for a = 1:nOmega

    Omega = Omega_list(a);

    root_folder = fullfile(base_folder, ...
        sprintf('Fresnel_Sweep_Omega=%.3f', Omega));

    for i = 1:nNF

        N_F = NF_list(i);

        folder_name = sprintf('NF_%0.1f', N_F);
        folder_name = strrep(folder_name, '.', 'p');

        workspace_file = fullfile( ...
            root_folder, folder_name, 'workspace.mat');

        if ~isfile(workspace_file)
            warning('Missing workspace: %s', workspace_file);
            continue;
        end

        data = load(workspace_file, 'sim', 'outputs');

        if isempty(sim_ref)
            sim_ref = data.sim;
        end

        % Common intensity normalization factor
        % if isfield(data.outputs, 'Gau_RHS')
        %     Gau_RHS{a,i} = abs(data.outputs.Gau_RHS).^2;
        % end
        % 
        % if isfield(data.outputs, 'Gau_LHS')
        %     Gau_LHS{a,i} = abs(data.outputs.Gau_LHS).^2;
        % end

    if isfield(data.outputs, 'Gau_RHS')
        I = abs(data.outputs.Gau_RHS).^2;
        if max(I, [], 'all') > 0
            Gau_RHS{a,i} = I / max(I, [], 'all');
        end
    end
    
    if isfield(data.outputs, 'Gau_LHS')
        I = abs(data.outputs.Gau_LHS).^2;
        if max(I, [], 'all') > 0
            Gau_LHS{a,i} = I / max(I, [], 'all');
        end
    end

    end
end

%% Plot settings

% Display window in meters
plot_limit = 0.6;

% Use one color scale per figure
cmap = parula(256);

%% SET 1: Fixed Fresnel number, varying Omega

for i = 1:nNF

    N_F = NF_list(i);

    fig = figure('Color','w', ...
        'Position',[100 100 2400 550], ...
        'Visible','off');

    tiledlayout(2, nOmega, ...
        'TileSpacing','compact', ...
        'Padding','compact');

    % Find maximum intensity across this figure
    Imax = 0;

    for a = 1:nOmega
        if ~isempty(Gau_RHS{a,i})
            Imax = max(Imax, max(Gau_RHS{a,i},[],'all'));
        end
        if ~isempty(Gau_LHS{a,i})
            Imax = max(Imax, max(Gau_LHS{a,i},[],'all'));
        end
    end

    if Imax == 0
        close(fig);
        continue;
    end

    for a = 1:nOmega

        % RHS mode
        nexttile(a);

        if ~isempty(Gau_RHS{a,i})
            imagesc(sim_ref.x, sim_ref.y, Gau_RHS{a,i});
            %clim([0 Imax]);
            clim([0 1]);
        end

        axis image;
        xlim([-plot_limit plot_limit]);
        ylim([-plot_limit plot_limit]);
        colormap(cmap);

        title(sprintf('\\Omega = %.3f rad/s', Omega_list(a)), ...
            'FontSize', 10);

        if a == 1
            ylabel('RHS   y [m]');
        end

        % LHS mode
        nexttile(nOmega+a);

        if ~isempty(Gau_LHS{a,i})
            imagesc(sim_ref.x, sim_ref.y, Gau_LHS{a,i});
            %clim([0 Imax]);
            clim([0 1]);
        end

        axis image;
        xlim([-plot_limit plot_limit]);
        ylim([-plot_limit plot_limit]);
        colormap(cmap);

        if a == 1
            ylabel('LHS   y [m]');
        end

        xlabel('x [m]');

    end

    sgtitle(sprintf( ...
        'Converged cavity modes | N_F = %.1f | Varying rotation rate', ...
        N_F), 'FontSize', 16, 'FontWeight', 'bold');

    filename = sprintf('Fixed_NF_%0.1f', N_F);
    filename = strrep(filename, '.', 'p');

    exportgraphics(fig, ...
        fullfile(save_folder, [filename '.png']), ...
        'Resolution', 300);

    savefig(fig, fullfile(save_folder, [filename '.fig']));

    close(fig);

    fprintf('Saved fixed NF = %.1f\n', N_F);

end

%% SET 2: Fixed Omega, varying Fresnel number

for a = 1:nOmega

    Omega = Omega_list(a);

    fig = figure('Color','w', ...
        'Position',[100 100 2400 550], ...
        'Visible','off');

    tiledlayout(2, nNF, ...
        'TileSpacing','compact', ...
        'Padding','compact');

    % Find maximum intensity across this figure
    Imax = 0;

    for i = 1:nNF
        if ~isempty(Gau_RHS{a,i})
            Imax = max(Imax, max(Gau_RHS{a,i},[],'all'));
        end
        if ~isempty(Gau_LHS{a,i})
            Imax = max(Imax, max(Gau_LHS{a,i},[],'all'));
        end
    end

    if Imax == 0
        close(fig);
        continue;
    end

    for i = 1:nNF

        % RHS mode
        nexttile(i);

        if ~isempty(Gau_RHS{a,i})
            imagesc(sim_ref.x, sim_ref.y, Gau_RHS{a,i});
            %clim([0 Imax]);
            clim([0 1]);
        end

        axis image;
        xlim([-plot_limit plot_limit]);
        ylim([-plot_limit plot_limit]);
        colormap(cmap);

        title(sprintf('N_F = %.1f', NF_list(i)), ...
            'FontSize', 10);

        if i == 1
            ylabel('RHS   y [m]');
        end

        % LHS mode
        nexttile(nNF+i);

        if ~isempty(Gau_LHS{a,i})
            imagesc(sim_ref.x, sim_ref.y, Gau_LHS{a,i});
            %clim([0 Imax]);
            clim([0 1]);
        end

        axis image;
        xlim([-plot_limit plot_limit]);
        ylim([-plot_limit plot_limit]);
        colormap(cmap);

        if i == 1
            ylabel('LHS   y [m]');
        end

        xlabel('x [m]');

    end

    sgtitle(sprintf( ...
        'Converged cavity modes | \\Omega = %.3f rad/s | Varying Fresnel number', ...
        Omega), 'FontSize', 16, 'FontWeight', 'bold');

    filename = sprintf('Fixed_Omega_%0.3f', Omega);
    filename = strrep(filename, '.', 'p');

    exportgraphics(fig, ...
        fullfile(save_folder, [filename '.png']), ...
        'Resolution', 300);

    savefig(fig, fullfile(save_folder, [filename '.fig']));

    close(fig);

    fprintf('Saved fixed Omega = %.3f\n', Omega);

end

fprintf('\nAll mode comparison figures saved.\n');