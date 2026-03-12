function [sim, laser, mirror, gain_medium] = Increase_Domain_Size(sim, laser, mirror, toggles, N_new)

% Increase the domain size for the first pass to mitigate boundary
% reflections.

[N_old, Ny] = size(laser.Gau);
if N_old ~= Ny
    error('Increase_Domain_Size: laser.Gau must be square.');
end
if N_new <= N_old
    error('Increase_Domain_Size: N_new must be > current field size.');
end

% Pad centered
G_new = zeros(N_new, N_new, 'like', laser.Gau);
c_old = floor(N_old/2) + 1;
c_new = floor(N_new/2) + 1;
idx_new = (c_new - (c_old-1)) : (c_new + (N_old - c_old));   % length N_old
G_new(idx_new, idx_new) = laser.Gau;
laser.Gau = G_new;

% Update sim struct
sim.N = N_new;
sim.Lx = sim.dx * sim.N;
x = (-floor(sim.N/2) : (ceil(sim.N/2)-1)) * sim.dx;
sim.x = x; sim.y = x;
[sim.X, sim.Y] = meshgrid(sim.x, sim.y);

% Absorbing mask
r = sqrt(sim.X.^2 + sim.Y.^2);
r0 = 0.8*(sim.Lx/2);
w = 0.1*(sim.Lx/2);
sim.mask_abs = ones(size(r));
ii = r > r0;
sim.mask_abs(ii) = exp(-((r(ii)-r0)/w).^8);

% Rebuild mirror masks on new grid
mirror(1).rmask = exp(-1i*laser.k0*(sim.X.^2+sim.Y.^2)/mirror(1).Rc);
mirror(2).rmask = exp(-1i*laser.k0*(sim.X.^2+sim.Y.^2)/mirror(2).Rc);
mirror(1).cmask = (sim.X.^2 + sim.Y.^2 <= (mirror(1).D/2)^2);
mirror(2).cmask = (sim.X.^2 + sim.Y.^2 <= (mirror(2).D/2)^2);

if toggles.gain_switch == true % resize gain medium to match
    gain_medium = Initialize_Gain_Medium(sim, mirror);
else
    gain_medium = [];
end

end