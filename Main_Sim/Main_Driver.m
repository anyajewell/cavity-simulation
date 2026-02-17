clear all;
clc;
figure;

% Simulation
[consts, sim, laser, frame, mirror, outputs, toggles] = Initialize_Sim(); % initialize
[gain_medium] = Initialize_Gain_Medium(sim, mirror);

%%
[laser, outputs, gain_medium] = Propagate_n_RTs(consts, sim, laser, frame, mirror, outputs, toggles, gain_medium); % propagate

%% Post-Processing
% Rescale 
% Pk(i) = max(max(abs(E).^2)); % peak intensity at this point
% E = E*sqrt(P0/sum(sum(abs(E).^2))); % rescale E
% BeamWidth(i) = sum(sum(sqrt(X.^2+Y.^2).*abs(E).^2))/P0;

% Analytic solution
x_rot = -frame.Omega/consts.c * (sim.z.^2 + 0.5*sim.L*sim.z);
x_acc = (frame.v0/consts.c)*(sim.z - sim.Z0) + (frame.accel/(2*consts.c^2))*(sim.z - sim.Z0).^2;
x_ana = x_rot + x_acc; % analytic solution from ray optics

% Plotting
%Plot_Results_vs_Analytic(centerx, x_ana, z) % plot results of propagation
Plot_Loss_Over_Prop(outputs.loss_frac, sim.RTs)
Plot_Center(outputs.centerx, sim.zs)

close(outputs.v);





