%% Load DeltaTheta study results from saved folders

clear; clc;

root_folder = "C:\Users\Anya Jewell\Documents\MATLAB\ORACLE\Results\2026-03-19\DeltaTheta_Loss_Study";

folders = dir(fullfile(root_folder, 'dtheta_x_*'));
folders = folders([folders.isdir]);

% Preallocate as cells
DeltaTheta_cell   = {};
lossRT_cell       = {};
nRT_used_cell     = {};
lossHist_cell     = {};
didConverge_cell  = {};
GauLHS_cell       = {};
GauRHS_cell       = {};
centerx_cell      = {};
centery_cell      = {};

k = 0;

for i = 1:numel(folders)
    folder_name = folders(i).name;

    % Parse DeltaTheta from folder name
    token = regexp(folder_name, 'dtheta_x_([-\d\.eE\+]+)', 'tokens', 'once');
    if isempty(token)
        warning('Could not parse DeltaTheta from folder: %s', folder_name);
        continue;
    end
    dtheta_val = str2double(token{1});

    ws_file = fullfile(root_folder, folder_name, 'workspace.mat');
    if ~isfile(ws_file)
        warning('Missing workspace.mat in folder: %s', folder_name);
        continue;
    end

    S = load(ws_file, 'outputs', 'sim', 'laser');
    sim = S.sim; laser = S.laser;

    k = k + 1;
    DeltaTheta_cell{k} = dtheta_val;

    % Converged RT loss
    if isfield(S.outputs, 'study') && isfield(S.outputs.study, 'loss_RT') && ~isempty(S.outputs.study.loss_RT)
        lossRT_cell{k} = S.outputs.study.loss_RT(end);
    elseif isfield(S.outputs, 'loss1') && isfield(S.outputs, 'loss2') ...
            && ~isempty(S.outputs.loss1) && ~isempty(S.outputs.loss2)
        lossRT_cell{k} = 1 - (1 - S.outputs.loss1(end))*(1 - S.outputs.loss2(end));
    else
        lossRT_cell{k} = NaN;
    end

    % Number of RTs used
    if isfield(S.outputs, 'study') && isfield(S.outputs.study, 'nRT_used') && ~isempty(S.outputs.study.nRT_used)
        nRT_used_cell{k} = S.outputs.study.nRT_used(end);
    elseif isfield(S.outputs, 'loss_frac')
        nRT_used_cell{k} = numel(S.outputs.loss_frac);
    else
        nRT_used_cell{k} = NaN;
    end

    % Optional convergence flag
    if isfield(S.outputs, 'did_converge')
        didConverge_cell{k} = S.outputs.did_converge;
    else
        didConverge_cell{k} = NaN;
    end

    % Loss history
    if isfield(S.outputs, 'loss_frac')
        lossHist_cell{k} = S.outputs.loss_frac(:);
    else
        lossHist_cell{k} = [];
    end

    % Converged mirror-plane fields
    if isfield(S.outputs, 'Gau_LHS')
        GauLHS_cell{k} = S.outputs.Gau_LHS;
    else
        GauLHS_cell{k} = [];
    end

    if isfield(S.outputs, 'Gau_RHS')
        GauRHS_cell{k} = S.outputs.Gau_RHS;
    else
        GauRHS_cell{k} = [];
    end

    % Centroid histories
    if isfield(S.outputs, 'centerx')
        centerx_cell{k} = S.outputs.centerx(:);
    else
        centerx_cell{k} = [];
    end

    if isfield(S.outputs, 'centery')
        centery_cell{k} = S.outputs.centery(:);
    else
        centery_cell{k} = [];
    end
end

% Convert scalar cells to arrays
DeltaTheta   = cell2mat(DeltaTheta_cell);
loss_RT      = cell2mat(lossRT_cell);
nRT_used     = cell2mat(nRT_used_cell);
did_converge = cell2mat(didConverge_cell);

% Sort by DeltaTheta
[DeltaTheta, idx] = sort(DeltaTheta);
loss_RT           = loss_RT(idx);
nRT_used          = nRT_used(idx);
did_converge      = did_converge(idx);

lossHist_cell     = lossHist_cell(idx);
GauLHS_cell       = GauLHS_cell(idx);
GauRHS_cell       = GauRHS_cell(idx);
centerx_cell      = centerx_cell(idx);
centery_cell      = centery_cell(idx);

%% Plot 1: converged loss vs DeltaTheta
figure;
mask = DeltaTheta > 0;
loglog(DeltaTheta(mask) / laser.theta_D, 100*loss_RT(mask), 'o-', 'LineWidth', 1.5);
hold on;
if any(DeltaTheta == 0) && any(mask)
    loglog(0, loss_RT(DeltaTheta==0), 'rs', 'MarkerFaceColor', 'r');
end
grid on;
xlabel('\Delta\theta_x / \theta_D');
ylabel('Converged loss / round-trip (%)');
title('Converged Loss vs Mirror Misalignment');

%% Plot 2: RTs used vs DeltaTheta
figure;
semilogx(DeltaTheta(mask), nRT_used(mask), 'o-', 'LineWidth', 1.5);
hold on;
if any(DeltaTheta == 0) && any(mask)
    plot(0, nRT_used(DeltaTheta==0), 'rs', 'MarkerFaceColor', 'r');
end
grid on;
xlabel('\Delta\theta_x [rad]');
ylabel('Round-trips used');
title('Reconvergence Time vs Mirror Misalignment');

%% Plot 3: transient RT loss histories
figure; hold on;
for i = 1:numel(DeltaTheta)
    if ~isempty(lossHist_cell{i})
        plot(1:numel(lossHist_cell{i}), lossHist_cell{i}, ...
            'DisplayName', sprintf('\\Delta\\theta = %.2e', DeltaTheta(i)), 'LineWidth', 1.5);
    end
end
grid on;
xlabel('Round-trip index');
ylabel('Loss / round-trip');
title('Transient Loss Histories');
legend('show', 'Location', 'eastoutside');

%% Plot 4: centroid x histories
figure; hold on;
for i = 1:numel(DeltaTheta)
    if ~isempty(centerx_cell{i})
        plot(1:numel(centerx_cell{i}), centerx_cell{i}, ...
            'DisplayName', sprintf('\\Delta\\theta = %.2e', DeltaTheta(i)));
    end
end
grid on;
xlabel('Sample index');
ylabel('Center x [m]');
title('Centroid x Histories');
legend('show', 'Location', 'eastoutside');

%% Plot 5: centroid y histories
figure; hold on;
for i = 1:numel(DeltaTheta)
    if ~isempty(centery_cell{i})
        plot(1:numel(centery_cell{i}), centery_cell{i}, ...
            'DisplayName', sprintf('\\Delta\\theta = %.2e', DeltaTheta(i)));
    end
end
grid on;
xlabel('Sample index');
ylabel('Center y [m]');
title('Centroid y Histories');
legend('show', 'Location', 'eastoutside');

%% Plot 6: converged Gau_LHS intensity for each case
for i = 1:numel(DeltaTheta)
    if ~isempty(GauLHS_cell{i})
        figure;
        imagesc(sim.X(1,:), sim.Y(:,1), abs(GauLHS_cell{i}).^2);
        axis image;
        %colorbar;
        title(sprintf('Converged LHS mode, \\Delta\\theta = %.2e rad', DeltaTheta(i)));
        axis([-0.5 0.5 -0.5 0.5]);
        xlabel('x [m]');
        ylabel('y [m]');
    end
end

%% Plot 7: converged Gau_RHS intensity for each case
for i = 1:numel(DeltaTheta)
    if ~isempty(GauRHS_cell{i})
        figure;
        imagesc(sim.X(1,:), sim.Y(:,1), abs(GauRHS_cell{i}).^2);
        axis image;
        %colorbar;
        title(sprintf('Converged RHS mode, \\Delta\\theta = %.2e rad', DeltaTheta(i)));
        axis([-0.5 0.5 -0.5 0.5]);
        xlabel('x [m]');
        ylabel('y [m]');
    end
end

%% Save compiled results
save(fullfile(root_folder, 'compiled_DeltaTheta_study.mat'), ...
    'DeltaTheta', 'loss_RT', 'nRT_used', 'did_converge', ...
    'lossHist_cell', 'GauLHS_cell', 'GauRHS_cell', ...
    'centerx_cell', 'centery_cell');

%% Save all open figures

save_fig_folder = fullfile(root_folder, 'compiled_figures');
if ~exist(save_fig_folder, 'dir')
    mkdir(save_fig_folder);
end

figHandles = findall(0, 'Type', 'figure');

for k = 1:length(figHandles)
    fig = figHandles(k);

    % Bring figure into focus (optional but safe)
    figure(fig);

    % Try to get a clean name
    fig_name = get(get(gca, 'Title'), 'String');
    
    if isempty(fig_name)
        fig_name = sprintf('figure_%d', k);
    end

    % Clean filename (remove spaces, special chars)
    fig_name = regexprep(fig_name, '[^\w\d]', '_');

    % Save as PNG
    saveas(fig, fullfile(save_fig_folder, [fig_name '.png']));

    % Also save MATLAB figure (optional but VERY useful)
    savefig(fig, fullfile(save_fig_folder, [fig_name '.fig']));
end