%% Initialization

clear all;
clc;
figure;

[consts, sim, laser, frame, mirror, outputs, toggles] = Initialize_Sim(); % initialize
[gain_medium] = Initialize_Gain_Medium(sim, mirror);

%% Propagation
if laser.pos ~= mirror(1).loc && laser.pos ~= mirror(2).loc % wavefront is starting within the cavity
    if sim.dz > 0
        [laser, outputs, sim, gain_medium] = R(consts, sim, laser, frame, mirror, outputs, toggles, gain_medium, 0);
    else
        [laser, outputs, sim, gain_medium] = L(consts, sim, laser, frame, mirror, outputs, toggles, gain_medium, 0);
    end
    [laser, outputs, gain_medium] = Propagate_n_RTs(consts, sim, laser, frame, mirror, outputs, toggles, gain_medium); % now propagate # of desired RTs
else % laser is starting at one end of the cavity
    [laser, outputs, gain_medium] = Propagate_n_RTs(consts, sim, laser, frame, mirror, outputs, toggles, gain_medium);
end

%% Post-Processing

% Analytic solution
% x_rot = -frame.Omega/consts.c * (outputs.zs.^2 + 0.5*sim.L*outputs.zs);
% x_acc = (frame.v0/consts.c)*(outputs.zs - sim.Z0) + (frame.accel/(2*consts.c^2))*(outputs.zs - sim.Z0).^2;
% x_ana = x_rot + x_acc; % analytic solution from ray optics
%Plot_Results_vs_Analytic(centerx, x_ana, z)

fig1 = Plot_RT_Loss(outputs.loss_frac, sim.RTs);
savefig(fig1, fullfile(outputs.saveFolder,'RT_Loss.fig'));
exportgraphics(fig1, fullfile(outputs.saveFolder,'RT_Loss.png'), 'Resolution',300);

fig2 = Plot_Center(outputs.centerx, outputs.zs);
savefig(fig2, fullfile(outputs.saveFolder,'Center_vs_z.fig'));
exportgraphics(fig2, fullfile(outputs.saveFolder,'Center_vs_z.png'), 'Resolution',300);

% fig3 = Plot_R(outputs.R1, outputs.R2, sim.Nz, outputs.zs);
% savefig(fig3, fullfile(outputs.saveFolder,'Reflection.fig'));
% exportgraphics(fig3, fullfile(outputs.saveFolder,'Reflection.png'), 'Resolution',300);
% 
% fig4 = Plot_Mirror_Loss(outputs.loss1, outputs.loss2, sim.Nz, outputs.zs);
% savefig(fig4, fullfile(outputs.saveFolder,'Mirror_Loss.fig'));
% exportgraphics(fig4, fullfile(outputs.saveFolder,'Mirror_Loss.png'), 'Resolution',300);

fig5 = Plot_Imax_vs_Isat(outputs, gain_medium.I_sat, sim);
savefig(fig5, fullfile(outputs.saveFolder,'Imax.fig'));
exportgraphics(fig5, fullfile(outputs.saveFolder,'Imax.png'), 'Resolution',300);

fig6 = Plot_RT_Gain(outputs.gain, sim.RTs);
savefig(fig6, fullfile(outputs.saveFolder,'Gain.fig'));
exportgraphics(fig6, fullfile(outputs.saveFolder,'Gain.png'), 'Resolution',300);

Save_Workspace(consts, sim, laser, frame, mirror, outputs, gain_medium); % automatic save

close(outputs.v);