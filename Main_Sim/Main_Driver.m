clear all;
clc;
figure;

% Simulation
[consts, sim, laser, frame, mirror, outputs, toggles] = Initialize_Sim(); % initialize
[gain_medium] = Initialize_Gain_Medium(sim, mirror);

%%
[laser, outputs, gain_medium] = Propagate_n_RTs(consts, sim, laser, frame, mirror, outputs, toggles, gain_medium); % propagate

%% Post-Processing

% Analytic solution
x_rot = -frame.Omega/consts.c * (sim.z.^2 + 0.5*sim.L*sim.z);
x_acc = (frame.v0/consts.c)*(sim.z - sim.Z0) + (frame.accel/(2*consts.c^2))*(sim.z - sim.Z0).^2;
x_ana = x_rot + x_acc; % analytic solution from ray optics

% Plotting
%Plot_Results_vs_Analytic(centerx, x_ana, z)
Plot_RT_Loss(outputs.loss_frac, sim.RTs)
Plot_Center(outputs.centerx, sim.zs)
Plot_R(outputs.R1, outputs.R2, sim.Nz, sim.zs)
Plot_Mirror_Loss(outputs.loss1, outputs.loss2, sim.Nz, sim.zs)

close(outputs.v);





