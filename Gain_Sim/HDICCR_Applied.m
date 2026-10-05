%% 3D steady-state thermal model: a-cut Nd:YVO4 rod, end-pumped
clear; close all

%% ---------------- Crystal geometry ----------------
R      = 1.0e-3;     % rod radius, m            [PLACEHOLDER] 2 mm dia carried over from previous model; get from spec sheet
L      = 5e-3;       % rod length, m            [PLACEHOLDER] typical length for ~1 at% doping; get from spec sheet
doping = 1.0;        % Nd concentration, at%    [PLACEHOLDER] common commercial doping; get from spec sheet
                     %   (not used in equations; alpha, k, eta_h below all ASSUME this doping)

%% ---------------- Thermal conductivity ----------------
k_source = "sato";   % "sato" (peer-reviewed measurement) or "datasheet" (vendor, lower bound)
switch k_source
    case "sato"
        kc = 12.1;   % along c (x), W/(m K)     [SOURCED/ESTIMATE] Sato & Taira [1], as quoted secondhand;
                     %   verify against [1] Table 2 and adjust for doping with their dk/dC_Nd (Fig. 9)
        ka = 8.9;    % along a (y,z), W/(m K)   [SOURCED/ESTIMATE] same caveat as kc
    case "datasheet"
        kc = 5.23;   % along c (x), W/(m K)     [SOURCED] vendor datasheets [2]; origin of value unclear
        ka = 5.10;   % along a (y,z), W/(m K)   [SOURCED] vendor datasheets [2]
end

isotropic_check = false;   % true -> kc = ka, for comparison against 2D axisymmetric model
if isotropic_check, kc = ka; end

% Temperature-dependent conductivity: k(T) = k_ref * (T_ref/T)^n_k, applied to all three axes
T_ref  = 298.15;     % reference temperature of kc, ka, K   [ESTIMATE] assumes kc, ka are room-temperature values
n_k    = 1;          % exponent                             [ESTIMATE] Umklapp-limited ~1/T; check Sato & Taira [1] Fig. 8
                     %   n_k = 0 -> constant k (must reproduce the previous model exactly)
compare_constant_k = false;  % true -> also solve with constant k and report the difference
run_convergence    = false;  % true -> run mesh convergence study at the end (slow; turn off once mesh is chosen)

%% ---------------- Other material properties ----------------
rho    = 4220;       % density, kg/m^3          [SOURCED] [2]; consistent across all sources (transient only)
cp     = 505;        % specific heat, J/(kg K)  [ESTIMATE] 24.6 cal/(mol K) from [3], converted with
                     %   M(YVO4) = 203.8 g/mol (transient only)
dne_dT = 2.9e-6;     % extraordinary dn/dT, 1/K [SOURCED] [3]; for later OPD step, pi-polarized beam

%% ---------------- Pump and heat deposition ----------------
lambda_p = 808e-9;   % pump wavelength, m       [PLACEHOLDER] standard diode wavelength; confirm with Anya
lambda_l = 1064.3e-9;% laser wavelength, m      [SOURCED] [3]
P      = 10;         % pump power at crystal, W [PLACEHOLDER] get measured power at crystal face
w      = 0.45e-3;    % pump 1/e^2 radius, m     [PLACEHOLDER] carried over from previous model, unconfirmed
alpha  = 2300;       % absorption coeff, 1/m    [SOURCED] 23 cm^-1 for 1.0 at% a-cut, pi-pol. peak [3];
                     %   likely an OVERESTIMATE for a broadband / unpolarized diode
eta_q  = 1 - lambda_p/lambda_l;   % quantum defect, ~0.24 at 808 nm (reference only)
eta_h  = 0.30;       % heat fraction            [ESTIMATE] non-lasing case. Lasing value 0.24 for 1 at% [4];
                     %   scaled by 1/0.77 from lens ratio [5] (assumes lens ~ heat load). Sensitivity range 0.24-0.40
T0     = 20;         % mount temperature, C     [PLACEHOLDER] get chiller/mount setpoint

Q0     = 2*eta_h*P*alpha/(pi*w^2);          % peak volumetric heat, W/m^3
P_heat = eta_h*P*(1 - exp(-alpha*L));       % total heat deposited, W

fprintf("k_source = %s: kc = %.2f, ka = %.2f W/(m K) | eta_h = %.2f | alpha = %.0f 1/m\n", ...
        k_source, kc, ka, eta_h, alpha);

%% Geometry
g = multicylinder(R, L);          % base at z = 0, top at z = L, axis along z

%% Model
model = femodel(AnalysisType="thermalSteady", Geometry=g);
model.MaterialProperties = materialProperties(ThermalConductivity=[kc ka ka]);   % constant k (replaced by k(T) before final solve)

%% Volumetric Gaussian pump source with Beer-Lambert absorption
model.CellLoad = cellLoad(Heat=@(loc,state) ...
    Q0*exp(-2*(loc.x.^2 + loc.y.^2)/w^2).*exp(-alpha*loc.z));

%% Boundary conditions
% F1 = z=0 pump face, F2 = z=L back face, F3 = side (verified)
% End faces left alone -> insulated
model.FaceBC(3) = faceBC(Temperature=T0);   % [ASSUMPTION] perfect thermal contact on full side

%% Initial guess (needed because the load is a function handle)
model.CellIC = cellIC(Temperature=T0);

%% Mesh (heat deposited within ~1/alpha, so keep elements small)
model = generateMesh(model, Hmax=min(w, 1/alpha)/3);
fprintf("Mesh: %d nodes, %d elements\n", size(model.Geometry.Mesh.Nodes,2), ...
        size(model.Geometry.Mesh.Elements,2));

%% Solve
% (a) Optional reference solve with constant k, on the same mesh
if compare_constant_k
    result_const = solve(model);
    Tpeak_const  = interpolateTemperature(result_const, [0;0;0]);   % same sampling as Tpeak
end

% (b) Final solve with temperature-dependent k.
%     state.u = local temperature (C); +273.15 converts to K for the ratio.
%     Returns a 3-by-N matrix: [kx; ky; kz] at every evaluation point.
kfun = @(loc,state) [kc; ka; ka] .* (T_ref ./ (state.u + 273.15)).^n_k;
model.MaterialProperties = materialProperties(ThermalConductivity=kfun);
result = solve(model);

%% Post-processing: sample the solution on grids
n = 201;
[X, Y] = meshgrid(linspace(-R, R, n));
Tface = reshape(interpolateTemperature(result, [X(:)'; Y(:)'; zeros(1,n^2)]), size(X));

zv = linspace(0, L, 301);  xv = linspace(-R, R, n);
[XZ, ZX] = meshgrid(xv, zv);
Txz = reshape(interpolateTemperature(result, [XZ(:)'; zeros(1,numel(XZ)); ZX(:)']), size(XZ));
Tyz = reshape(interpolateTemperature(result, [zeros(1,numel(XZ)); XZ(:)'; ZX(:)']), size(XZ));

r  = linspace(0, R, 200);
Tx = interpolateTemperature(result, [r; zeros(size(r)); zeros(size(r))]);
Ty = interpolateTemperature(result, [zeros(size(r)); r; zeros(size(r))]);

Qout  = evaluateHeatRate(result, "Face", 3);
Tpeak = interpolateTemperature(result, [0;0;0]);   % T at pump-face center (true peak location)
clim_all = [T0 Tpeak];            % shared color scale for all heat maps

%% All plots in one figure (2 x 3 grid)
fig = figure(Name="Nd:YVO4 a-cut thermal model", Units="normalized", Position=[0.03 0.08 0.94 0.82]);
tl = tiledlayout(fig, 2, 3, TileSpacing="compact", Padding="compact");

% (1) Geometry with face labels
nexttile(tl, 1)
pdegplot(g, FaceLabels="on", FaceAlpha=0.4)
title("Geometry (F1 pump face, F3 side = mount)")

% (2) Pump face, face-on
nexttile(tl, 2)
imagesc(X(1,:)*1e3, Y(:,1)*1e3, Tface, AlphaData=~isnan(Tface))
axis xy equal tight; clim(clim_all); colorbar
hold on; contour(X*1e3, Y*1e3, Tface, 10, "k"); hold off
xlabel("x, c-axis (mm)"); ylabel("y, a-axis (mm)")
title("Pump face (z = 0), T (C)")

% (3) Radial profiles along c and a
nexttile(tl, 3)
plot(r*1e3, Tx - T0, r*1e3, Ty - T0, "--")
xlabel("Distance from axis (mm)"); ylabel("\DeltaT (K)")
legend("Along c (x)", "Along a (y)"); grid on
title("Pump-face profile: c vs a")

% (4) Longitudinal slice along c
nexttile(tl, 4)
imagesc(zv*1e3, xv*1e3, Txz'); axis xy tight; clim(clim_all); colorbar
xlabel("z (mm)"); ylabel("x, c-axis (mm)")
title("Longitudinal slice along c (y = 0)")

% (5) Longitudinal slice along a
nexttile(tl, 5)
imagesc(zv*1e3, xv*1e3, Tyz'); axis xy tight; clim(clim_all); colorbar
xlabel("z (mm)"); ylabel("y, a-axis (mm)")
title("Longitudinal slice along a (x = 0)")

% (6) Anisotropy: difference between c and a profiles
nexttile(tl, 6)
plot(r*1e3, Tx - Ty); grid on
xlabel("Distance from axis (mm)"); ylabel("T_c - T_a (K)")
title("Anisotropy: T_c - T_a at pump face")

% Summary line across the top
title(tl, sprintf(['k: %s (k_c = %.2f, k_a = %.2f W/m K at T_{ref}), n_k = %g  |  \\eta_h = %.2f  |  ' ...
                   'Peak \\DeltaT = %.1f K  |  Heat in %.3f W / out %.3f W'], ...
                   k_source, kc, ka, n_k, eta_h, Tpeak - T0, P_heat, abs(Qout)));

%% Summary to command window
fprintf("Peak T: %.2f C (dT = %.2f K) with k(T), n_k = %g\n", Tpeak, Tpeak - T0, n_k);
if compare_constant_k
    fprintf("Peak T with constant k: %.2f C (dT = %.2f K) | k(T) adds %.2f K (%+.1f%%)\n", ...
            Tpeak_const, Tpeak_const - T0, Tpeak - Tpeak_const, 100*(Tpeak - Tpeak_const)/(Tpeak_const - T0));
end
fprintf("Heat deposited: %.3f W | Heat out through side: %.3f W\n", P_heat, abs(Qout));

%% Mesh convergence study
% Re-meshes the full model (k(T), source, BCs included) at increasing resolution
% and tracks dT at the pump-face center. Converged when the last two agree within ~1%.
if run_convergence
    divs = [3 4 5];                        % Hmax = min(w,1/alpha)/divs
    Tpk = zeros(size(divs)); Nn = Tpk;
    for i = 1:numel(divs)
        m   = generateMesh(model, Hmax=min(w,1/alpha)/divs(i));
        res = solve(m);
        Tpk(i) = interpolateTemperature(res, [0;0;0]) - T0;   % dT at pump-face center
        Nn(i)  = size(m.Geometry.Mesh.Nodes, 2);
        fprintf("Hmax = 1/%d: %7d nodes, peak dT = %.3f K\n", divs(i), Nn(i), Tpk(i));
    end
    fprintf("Change between last two meshes: %.2f%%\n", 100*abs(Tpk(end) - Tpk(end-1))/Tpk(end));

    figure(Name="Mesh convergence")
    plot(Nn, Tpk, "o-"); grid on
    xlabel("Mesh nodes"); ylabel("Peak \DeltaT (K)"); title("Mesh convergence")
end