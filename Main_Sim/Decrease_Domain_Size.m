function [sim, laser] = Decrease_Domain_Size(sim, laser)

% Restore domain to original initialization size.

if ~isfield(sim, 'orig')
    error('Original domain parameters not stored.');
end

N_target = sim.orig.N;
N_current = sim.N;

if N_target >= N_current
    error('Target domain must be smaller than current.');
end

center = floor(N_current/2) + 1;
half = floor(N_target/2);
idx = (center-half):(center+half-1);
laser.Gau = laser.Gau(idx, idx);

% Restore sim struct
sim.N = sim.orig.N;
sim.Lx = sim.orig.Lx;
sim.dx = sim.orig.dx;
x = (-sim.N/2 : sim.N/2-1) * sim.dx;
sim.x = x; sim.y = x;
[sim.X, sim.Y] = meshgrid(sim.x, sim.y);

end