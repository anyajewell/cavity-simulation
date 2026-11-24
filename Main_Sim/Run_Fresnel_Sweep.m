clear;
clc;
close all;

NF_list = 0.1:0.1:1.0;
Omega_list = 0;
total_cases = numel(NF_list) * numel(Omega_list);

base_folder = 'C:\Users\jewellan\Documents\GitHub\cavity-simulation\Results\SciTech_Data';

for a = 1:length(Omega_list)

    Omega = Omega_list(a);

    root_folder = fullfile(base_folder, ...
    sprintf('Fresnel_Sweep_Omega=%.3f', Omega));

    if ~exist(root_folder, 'dir')
        mkdir(root_folder);
    end

    for i = 1:length(NF_list)
        close all;

        case_num = (a-1)*numel(NF_list) + i;
    
        N_F = NF_list(i);

        fprintf('\n====================================\n');
        fprintf('Running case %d/%d (%.1f%%)\n', ...
            case_num, total_cases, 100*case_num/total_cases);
        fprintf('Omega = %.4f rad/s, NF = %.2f\n', Omega, N_F);
        fprintf('====================================\n');
    
        % Create unique folder for this Fresnel number
        folder_name = sprintf('NF_%0.1f', N_F);
        folder_name = strrep(folder_name, '.', 'p');
        case_folder = fullfile(root_folder, folder_name);
    
        if ~exist(case_folder, 'dir')
            mkdir(case_folder);
        end
    
        % Initialize everything through your function
        [consts, sim, laser, frame, mirror, outputs, toggles, gain_medium] = Initialize_Sim(N_F, Omega);
    
        % Redirect saving to this case folder
        outputs.saveFolder = case_folder;
    
        % Run simulation
        [laser, outputs, gain_medium] = Propagate_n_RTs( ...
            consts, sim, laser, frame, mirror, outputs, toggles, gain_medium);
    
        % Save workspace
        save(fullfile(case_folder, 'workspace.mat'), ...
            'consts', 'sim', 'laser', 'frame', 'mirror', 'outputs', ...
            'toggles', 'gain_medium', 'N_F');
    
        fprintf('Saved case to:\n%s\n', case_folder);
    end

end

fprintf('\nSweep complete.\n');