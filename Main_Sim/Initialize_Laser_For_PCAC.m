function [consts, sim, laser, frame, mirror, outputs, toggles, gain_medium, loss_frac] = Initialize_Laser_For_PCAC(path, file)

    fullpath = fullfile(path, file);
    [consts, sim, laser, frame, mirror, outputs, gain_medium, toggles] = Load_Workspace(fullpath); % load pre-solved cavity mode solution
    toggles.finish_line = 'convergence'; toggles.track_center = 0; toggles.outputs_switch = 0; % ensure no unnecessary outputs
    loss_frac = Get_Latest_Loss_Fraction(fullpath);

end
    