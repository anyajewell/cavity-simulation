clear; close all
load("Copper_results.mat")   % xg, zs, T3, T0, dne_dT, R, L, lambda_L, lasing, P, P_out, P_heat, kc, ka, n_k

%% ---------------- Inputs ----------------
n_pass  = 1;          % passes of the probe through the crystal   [PLACEHOLDER] confirm with Anya
a_fit   = 0.75e-3;    % fit aperture radius at the crystal, m     [CHOSEN] = pump waist w0; must match the WFS analysis
M_relay = 2.14;       % relay magnification, crystal to WFS31     [PLACEHOLDER]
lambda  = lambda_L;   % probe wavelength, m (1064 nm, pi-polarized, sees n_e)
if a_fit >= R, error("a_fit must be smaller than the crystal radius R = %.2f mm", R*1e3); end

%% ---------------- (1) Optical path difference ----------------
OPD = n_pass*trapz(zs, dne_dT*(T3 - T0), 3);     % m, on the (x, y) grid; NaN outside the crystal

%% ---------------- (2) Zernike fit over the aperture ----------------
[X, Y] = ndgrid(xg, xg);
r  = hypot(X, Y);  th = atan2(Y, X);
in = r <= a_fit & ~isnan(OPD);
Z  = zernike_osa(r(in)/a_fit, th(in));           % one column per j = 0..14
c  = Z \ OPD(in);                                % least-squares coefficients, m
rms_res = sqrt(mean((OPD(in) - Z*c).^2));        % fit residual, m
cw = c/lambda;                                   % coefficients in waves

%% ---------------- (3)-(5) Thermal lens ----------------
c4 = c(5);  c5 = c(6);                           % MATLAB index = j + 1
f_c = -a_fit^2/(2*(2*sqrt(3)*c4 + sqrt(6)*c5));  % along c (x)
f_a = -a_fit^2/(2*(2*sqrt(3)*c4 - sqrt(6)*c5));  % along a (y)

%% ---------------- Printout ----------------
if lasing, mode_str = sprintf("lasing, P_out = %.2f W", P_out); else, mode_str = "non-lasing"; end
fprintf("Input: P = %.1f W, %s, heat %.3f W | dn_e/dT = %.2e /K | %d pass(es), fit radius %.3f mm (%d points)\n", ...
        P, mode_str, P_heat, dne_dT, n_pass, a_fit*1e3, nnz(in));
fprintf("Peak OPD = %.3f waves (%.1f nm) at 1064 nm\n", max(OPD(:))/lambda, max(OPD(:))*1e9);
fprintf("Zernike (OSA, waves): Z3 %+.4f  Z4 %+.4f  Z5 %+.4f  Z7 %+.4f  Z8 %+.4f  Z12 %+.4f | residual rms %.4f\n", ...
        cw(4), cw(5), cw(6), cw(8), cw(9), cw(13), rms_res/lambda);
fprintf("Thermal lens at the crystal: f_c = %.3f m, f_a = %.3f m, f_c/f_a = %.3f\n", f_c, f_a, f_c/f_a);
fprintf("At the WFS31 (M = %.2f): f_c = %.2f m, f_a = %.2f m\n", M_relay, M_relay^2*f_c, M_relay^2*f_a);

%% ---------------- Plots ----------------
W_ap = nan(size(OPD));                           % wavefront in the aperture, piston and tilt removed
W_ap(in) = (OPD(in) - Z(:,1:3)*c(1:3))/lambda;
ph = linspace(0, 2*pi, 200);

figure(Name="Wavefront", Units="normalized", Position=[0.05 0.2 0.9 0.45]);
tiledlayout(1, 3, TileSpacing="compact", Padding="compact");
nexttile
imagesc(xg*1e3, xg*1e3, OPD'/lambda, AlphaData=~isnan(OPD')); axis xy equal tight; colorbar
hold on; plot(a_fit*1e3*cos(ph), a_fit*1e3*sin(ph), "w"); hold off
xlabel("x, c-axis (mm)"); ylabel("y, a-axis (mm)"); title("OPD over the crystal (waves)")
nexttile
imagesc(xg*1e3, xg*1e3, W_ap', AlphaData=~isnan(W_ap')); axis xy equal; colorbar
xlim([-1 1]*a_fit*1e3); ylim([-1 1]*a_fit*1e3)
xlabel("x, c-axis (mm)"); ylabel("y, a-axis (mm)"); title("Wavefront in aperture, piston/tilt removed (waves)")
nexttile
bar(3:14, cw(4:15)); grid on; xticks(3:14)
xlabel("OSA index j"); ylabel("Coefficient (waves)")
title(sprintf("Zernike terms | f_c = %.2f m, f_a = %.2f m", f_c, f_a))

%% ================= Local function =================
function Z = zernike_osa(rho, t)
% OSA/ANSI Zernike polynomials Z_j, j = 0..14, as listed in [Th] Table 2.
% rho = normalized radius (0..1), t = angle from +x (rad); column vectors. Returns N x 15.
Z = [ones(size(rho)), ...                                                       % j = 0  piston
     2*rho.*sin(t), 2*rho.*cos(t), ...                                          % 1, 2   tilt
     sqrt(6)*rho.^2.*sin(2*t), sqrt(3)*(2*rho.^2 - 1), sqrt(6)*rho.^2.*cos(2*t), ...   % 3, 4, 5
     sqrt(8)*rho.^3.*sin(3*t), sqrt(8)*(3*rho.^3 - 2*rho).*sin(t), ...          % 6, 7
     sqrt(8)*(3*rho.^3 - 2*rho).*cos(t), sqrt(8)*rho.^3.*cos(3*t), ...          % 8, 9
     sqrt(10)*rho.^4.*sin(4*t), sqrt(10)*(4*rho.^4 - 3*rho.^2).*sin(2*t), ...   % 10, 11
     sqrt(5)*(6*rho.^4 - 6*rho.^2 + 1), ...                                     % 12     spherical
     sqrt(10)*(4*rho.^4 - 3*rho.^2).*cos(2*t), sqrt(10)*rho.^4.*cos(4*t)];      % 13, 14
end
