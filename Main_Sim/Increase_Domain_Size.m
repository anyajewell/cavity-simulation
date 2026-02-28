function [sim, laser] = Increase_Domain_Size(sim, laser, N_new, mirror)
% Increase domain size by zero-padding laser.Gau to N_new x N_new.
% Robust to sim.N not matching size(laser.Gau). Keeps dx constant.

% --- Determine current field size from laser.Gau (source of truth) ---
[Nx, Ny] = size(laser.Gau);
if Nx ~= Ny
    error('Increase_Domain_Size: laser.Gau must be square, got %d-by-%d.', Nx, Ny);
end
N_old = Nx;

if N_new <= N_old
    error('Increase_Domain_Size: N_new (%d) must be > current field size (%d).', N_new, N_old);
end

% --- Save original domain once (from current sim and/or field) ---
if ~isfield(sim, 'orig') || ~isfield(sim.orig, 'N')
    sim.orig.Lx = sim.Lx;
    sim.orig.N  = N_old;     % store what the field actually is
    sim.orig.dx = sim.dx;
end

% --- Allocate new field ---
G_new = zeros(N_new, N_new, 'like', laser.Gau);

% --- Center-align using index ranges that always match N_old ---
c_old = floor(N_old/2) + 1;
c_new = floor(N_new/2) + 1;

i_new = (c_new - (c_old-1)) : (c_new + (N_old - c_old));  % length N_old

if numel(i_new) ~= N_old
    error('Increase_Domain_Size: internal index mismatch (expected %d, got %d).', N_old, numel(i_new));
end

G_new(i_new, i_new) = laser.Gau;

% --- Commit ---
laser.Gau = G_new;

% --- Update sim to match new field ---
sim.N  = N_new;
sim.Lx = sim.dx * sim.N;

x = (-floor(sim.N/2) : (ceil(sim.N/2)-1)) * sim.dx;  % works for odd/even
sim.x = x;
sim.y = x;
[sim.X, sim.Y] = meshgrid(sim.x, sim.y);

r = sqrt(sim.X.^2 + sim.Y.^2);
r0 = 0.8*(sim.Lx/2); % damping starts at 80% radius
w = 0.1*(sim.Lx/2); % damping width
sim.mask_abs = ones(size(r));
idx = r > r0;
sim.mask_abs(idx) = exp(-((r(idx)-r0)/w).^8); % steep super-Gaussian

% Mirror masks: reflecting lens phase screens and clipping masks
rmask1 = exp(-1i*laser.k0*(sim.X.^2+sim.Y.^2)/(mirror(1).Rc)); % reflection mask mirror 1 (RHS), negative sign because propagation uses positive sign convention
rmask2 = exp(-1i*laser.k0*(sim.X.^2+sim.Y.^2)/(mirror(2).Rc)); % reflection mask mirror 2 (LHS), negative sign because propagation uses positive sign convention
cmask1 = (sim.X.^2 + sim.Y.^2 <= (mirror(1).D/2)^2); % clipping mask mirror 1 (RHS)
cmask2 = (sim.X.^2 + sim.Y.^2 <= (mirror(2).D/2)^2); % clipping mask mirror 2 (LHS)
mirror(1).rmask = rmask1; mirror(2).rmask = rmask2; mirror(1).cmask = cmask1; mirror(2).cmask = cmask2;

end