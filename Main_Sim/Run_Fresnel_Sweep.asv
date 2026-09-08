clear;
clc;
close all;

NF_list = 0.1:0.1:1.0;

% Root folder for the whole sweep
root_folder = fullfile(pwd, 'Fresnel_Sweep');
if ~exist(root_folder, 'dir')
    mkdir(root_folder);
end

for i = 1:length(NF_list)
    close all;

    N_F = NF_list(i);

    fprintf('\n=============================\n');
    fprintf('Running case %d of %d: N_F = %.1f\n', i, length(NF_list), N_F);
    fprintf('=============================\n');

    % Create unique folder for this Fresnel number
    folder_name = sprintf('NF_%0.1f', N_F);
    folder_name = strrep(folder_name, '.', 'p');
    case_folder = fullfile(root_folder, folder_name);

    if ~exist(case_folder, 'dir')
        mkdir(case_folder);
    end

    % Initialize everything through your function
    [consts, sim, laser, frame, mirror, outputs, toggles, gain_medium] = Initialize_Sim(N_F);

    % Redirect saving to this case folder
    outputs.saveFolder = case_folder;

    % If video is on, reopen the video writer in this folder so files do not overwrite
    if toggles.outputs_switch == true
        try
            % Close the one created inside Initialize_Sim, if it exists
            if isfield(outputs, 'v')
                close(outputs.v);
            end
        catch
        end

        videoname = sprintf('Omega=%.3f_accel=%.2e_L=%.0f_D=%.2f_RTs=%.0f_NF=%.1f.mp4', ...
            frame.Omega, frame.accel, sim.L, mirror(1).D, sim.RTs, N_F);
        videoname = strrep(videoname, '.', 'p');

        outputs.v = VideoWriter(fullfile(case_folder, videoname), 'MPEG-4');
        open(outputs.v);
    end

    % Run simulation
    [laser, outputs, gain_medium] = Propagate_n_RTs( ...
        consts, sim, laser, frame, mirror, outputs, toggles, gain_medium);

    % Save workspace
    save(fullfile(case_folder, 'workspace.mat'), ...
        'consts', 'sim', 'laser', 'frame', 'mirror', 'outputs', ...
        'toggles', 'gain_medium', 'N_F');

    % Save open figures into this folder
    figs = findall(0, 'Type', 'figure');
    for k = 1:length(figs)
        fig_png = fullfile(case_folder, sprintf('fig_%02d.png', k));
        fig_fig = fullfile(case_folder, sprintf('fig_%02d.fig', k));

        try
            exportgraphics(figs(k), fig_png, 'Resolution', 300);
        catch
            saveas(figs(k), fig_png);
        end

        try
            savefig(figs(k), fig_fig);
        catch
        end
    end

    % Close video file for this case
    if toggles.outputs_switch == true
        try
            close(outputs.v);
        catch
        end
    end

    fprintf('Saved case to:\n%s\n', case_folder);
end

fprintf('\nFresnel sweep complete.\n');