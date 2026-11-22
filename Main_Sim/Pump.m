function [gain, F, N2, I, IL] = Pump(P_circ, p)

% Propagates the pump through the crystal with the laser field present and returns the round-trip gain of the laser mode, (5).
% Along each ray (fixed rho), F obeys (3): dF/dz = -sigma*(N - N2)*F. The intensity is
% I = (w(0)/w(z))^2 * F, so the power between two rays changes only by absorption (DERIVED from the
% 1/w^2 peak scaling of a Gaussian beam). F(rho, 0) = 2P/(pi w(0)^2) exp(-2 rho^2) carries power P.

    rho = p.rho;  zg = p.zg;  dz = zg(2) - zg(1);
    w_0 = p.wz(0);
    IL_at = @(z) 2*P_circ*p.u_L(rho*p.wz(z)); % forward + backward laser intensity
    dFdz = @(z, Fk) -p.sigma*(p.N - p.N2fun((w_0/p.wz(z))^2*Fk, IL_at(z))).*Fk;
    F = zeros(numel(rho), numel(zg));
    F(:,1) = 2*p.P/(pi*w_0^2)*exp(-2*rho.^2);
    
    for k = 1:numel(zg) - 1 % RK4 in z
        z = zg(k);  Fk = F(:,k);
        k1 = dFdz(z, Fk);
        k2 = dFdz(z + dz/2, Fk + dz/2*k1);
        k3 = dFdz(z + dz/2, Fk + dz/2*k2);
        k4 = dFdz(z + dz,   Fk + dz*k3);
        F(:,k+1) = Fk + dz/6*(k1 + 2*k2 + 2*k3 + k4);
    end
    
    wzg = p.wz(zg);
    I = (w_0./wzg).^2.*F;
    IL = 2*P_circ*p.u_L(rho.*wzg);
    N2 = p.N2fun(I, IL);
    gain = 2*trapz(zg, trapz(rho, p.sigma_L*N2.*p.u_L(rho.*wzg).*(2*pi*wzg.^2.*rho), 1));

end