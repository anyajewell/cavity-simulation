%% NdYVO4_dual_rod_CW_model.m

clear; clc; close all;

outdir = fullfile(fileparts(mfilename('fullpath')), 'figures');
if ~exist(outdir,'dir'); mkdir(outdir); end

%% ===================== 1. PHYSICAL CONSTANTS =====================
h = 6.62607015e-34;     % Planck constant [J s]
c = 2.99792458e8;       % speed of light [m/s]

%% ===================== 2. LASER / CRYSTAL PARAMETERS (Table 1, Shen et al. 2017) =====================

% --- wavelengths ---
lambda_p = 808e-9;      % pump wavelength [m]. Paper uses 888e-9 (in-band);
                        % switched here to 808 nm, the standard diode-peak
                        % pumping wavelength for Nd:YVO4. See Section 6.
lambda_s = 1064e-9;     % laser wavelength [m]                (Table 1: lambda_s)
lambda_f = 1032e-9;     % mean fluorescence wavelength [m]    (Table 1: lambda_f, ref [6])
nu_p = c/lambda_p;
nu_s = c/lambda_s;

% --- thermal properties ---
T0_amb = 298;           % ambient / heat-sink temperature [K], 25 C    (Table 1: T0)
Ka  = 5.10;             % thermal conductivity, a-axis [W/(m K)]       (Table 1: Ka, ref [16])
Kc  = 5.23;             % thermal conductivity, c-axis [W/(m K)]       (Table 1: Kc, ref [16])
Hconv = 2e4;            % side-face convective coefficient [W/(m^2 K)] (Table 1: H, ref [27])

% --- cavity / crystal geometry ---
l   = 0.63;             % full cavity length [m]                       (Table 1: l)
lm  = 0.06;             % AOM length [m]                               (Table 1: lm)
lc  = 0.025;            % length of EACH crystal rod [m]               (Table 1: lc)
w_x = 2e-3;             % crystal width [m]                            (Table 1: w)
hgt = 2e-3;             % crystal height [m]                           (Table 1: h)

% --- spectroscopy: ion density, lifetimes, ETU ---
n_tot = 6.24e25;        % total Nd3+ ion density [1/m^3], 0.5 at.% doping (Table 1: ntot)
Wup   = 0.8e-21;        % ETU (energy transfer upconversion) coefficient [m^3/s] (Table 1: Wup, ref [28])

tau1  = 530e-12;        % lifetime of level 1, 4I11/2 -- the LOWER laser level [s] (Table 1: tau1, ref [29])
tau4  = 104.29e-6;      % lifetime of level 4, 4F3/2 -- the UPPER laser level [s]  (Table 1: tau4, ref [10])
tauup = 20e-9;          % relaxation lifetime of the transient ETU-manifold population [s] (Table 1: tauup, ref [22])

% Fluorescence branching ratios:

beta40 = 0.420; beta41 = 0.467; beta42 = 0.110; beta43 = 0.003;   % (Table 1, ref [30])
beta4to1 = beta41 + beta42 + beta43;

% --- cross sections ---
sigma_ap0 = 2.484e-23;  % PUMP ABSORPTION cross section [m^2] @ 808 nm
                        
                        
                        
sigma_es0 = 15.6e-23;   % LASER EMISSION cross section [m^2] @ 1064 nm  (Table 1: sigma_es0, ref [31]).
                       

% --- refractive indices ---
no = 1.96; ne = 2.17;   % ordinary/extraordinary refractive indices @ 1064nm
nc_ref = (no+ne)/2;     % mean index used for the a-cut crystal approximation (paper, Sec. 2)
nm = 1.45;              % AOM refractive index                          (Table 1: nm, ref [25])

% --- cavity losses / mode sizes ---
Toc  = 0.57;            % output coupler transmission          (Table 1: TOC)
Lloss = 0.02;           % intrinsic (distributed) cavity loss  (Table 1: L)

omega0   = 420e-6;      % laser mode radius [m]                (Table 1: omega0)
omega_p0 = 450e-6;      % pump mode radius [m]                 (Table 1: omega_p0)

%% ===================== Derived cavity quantities =====================
Vm   = pi*omega0^2 * l;

leff = l + 2*lc*(nc_ref-1) + lm*(nm-1);
tau_r = 2*leff/c;

fprintf('Vm = %.4g m^3, leff = %.4f m, tau_r = %.4g s (%.2f ns)\n', Vm, leff, tau_r, tau_r*1e9);

%% ===================== 3. PUMP RATE (Eq. 8) =====================
V_gain_each  = pi*omega_p0^2 * lc;                % pumped volume of ONE crystal [m^3]
alpha0       = sigma_ap0 * n_tot;                 % absorption coefficient [1/m] = sigma_ap*n_tot
absfrac      = 1 - exp(-alpha0*lc);               % fraction of pump absorbed over one crystal (Beer's law)
N_ions_each  = n_tot * V_gain_each;                % number of active ions in the pumped volume of one crystal

Rp_of = @(Pp_each) (Pp_each*absfrac) / (h*nu_p*N_ions_each);   % per-ion pump rate [1/s], Eq. 8

fprintf('alpha0 = %.1f 1/m, absorption depth 1/alpha0 = %.3f mm, absfrac = %.4f\n', ...
        alpha0, 1/alpha0*1e3, absfrac);

%% ===================== 4. CW STEADY STATE =====================
Dn_th = (1/tau_r)*(Lloss + log(1/(1-Toc))) / (2*(lc/leff)*(c/nc_ref)*sigma_es0);   % Eq. 14, solved for threshold
fprintf('Threshold inversion Dn_th = %.4g m^-3\n', Dn_th);

Pp_each_CW = 67;     % W per crystal (134 W total pump, the paper's CW test point, Sec. 4.1)
[Ps_cw, phi_cw, n4_cw, n1_cw] = cw_state(Pp_each_CW, Rp_of, Dn_th, n_tot, tau4, tau1, ...
                                          beta4to1, Wup, nc_ref, c, sigma_es0, ...
                                          nu_s, Vm, tau_r, Toc);
Ws_cw = (c/nc_ref)*sigma_es0*phi_cw;   % Eq. 9: Ws = (c/nc)*sigma_es*phi_s

fprintf('\n--- CW steady state, 134 W total pump ---\n');
fprintf('n4 = %.4g m^-3, n1 = %.4g m^-3, Dn = %.4g m^-3\n', n4_cw, n1_cw, n4_cw-n1_cw);
fprintf('Simulated output power = %.2f W   (paper: 62.8 W sim / 61.6 W expt AT 888NM -- not directly comparable, see header note)\n', Ps_cw);

%% ===================== 5. HEAT GENERATION (Eq. 18) =====================
E1_avg_cm1 = mean([1966 1988 2047 2062 2154 2182]);   % 4I11/2 sublevels, Y1-Y6 (Fig. 1)
E0_avg_cm1 = mean([433 226 173 108 0]);               % 4I9/2 sublevels, Z1-Z5 (Fig. 1)
E4_avg_cm1 = mean([11366 11384]);                     % 4F3/2 sublevels, R1,R2 (Fig. 1)

E10 = wn2J(E1_avg_cm1 - E0_avg_cm1, h, c);              % heat per n1->n0 decay
E41 = wn2J(E4_avg_cm1 - E1_avg_cm1, h, c);              % heat per ETU-manifold relaxation event
Ef  = h*c/lambda_f;
E4f = wn2J(E4_avg_cm1, h, c) - Ef;                       % quantum-defect heat per n4 decay

% ---------------------------------------------------------------------
% THE 888nm -> 808nm TRANSITION, IN FULL:
E_R1 = wn2J(11366, h, c);           % R1 sublevel energy
E_pump_defect = max(h*c/lambda_p - E_R1, 0);   % extra heat per absorbed photon; ~0 for 888nm (near-resonant), >0 for 808nm
fprintf('Extra pump-relaxation heat/absorbed photon = %.3e J (%.1f%% of pump photon energy)\n', ...
        E_pump_defect, E_pump_defect/(h*c/lambda_p)*100);
% ---------------------------------------------------------------------

V_xtal = w_x*hgt*lc;
Aside  = 2*(w_x+hgt)*lc;   % four convectively-cooled side faces (end faces adiabatic, Eqs. 21-25)

%% ===================== 6. PLOTS 1 & 2: CW POWER & TEMPERATURE VS. PUMP POWER =====================
Pp_each_list = linspace(15,70,25);
Ps_list = zeros(size(Pp_each_list));
T_list  = zeros(size(Pp_each_list));
for k = 1:numel(Pp_each_list)
    [Ps_list(k), phik, n4k, n1k] = cw_state(Pp_each_list(k), Rp_of, Dn_th, n_tot, tau4, tau1, ...
                                             beta4to1, Wup, nc_ref, c, sigma_es0, ...
                                             nu_s, Vm, tau_r, Toc);
    nupk = Wup*n4k^2*tauup;                                    % Eq. 1, adiabatic quasi-steady value
    Qk = E10*n1k/tau1 + E41*nupk/tauup + E4f*n4k/tau4 ...
         + E_pump_defect*Rp_of(Pp_each_list(k))*(n_tot-n1k-n4k);   % Eq. 18 + 808nm pump-relaxation term
    T_list(k) = (T0_amb-273.15) + Qk*V_xtal/(Hconv*Aside);     % Eq. 19 + Eqs. 21-24, lumped
end

save_plot(outdir, 'cw_power_vs_pump', 2*Pp_each_list, Ps_list, ...
          'Total pump power [W]', 'Output power [W]', ...
          'CW output power vs. pump power (cf. paper Fig. 8a)');

save_plot(outdir, 'cw_temperature_vs_pump', 2*Pp_each_list, T_list, ...
          'Total pump power [W]', 'Average crystal temperature [C]', ...
          'CW average crystal temperature vs. pump power (cf. paper Fig. 8b)');

Qcw = E10*n1_cw/tau1 + E41*(Wup*n4_cw^2*tauup)/tauup + E4f*n4_cw/tau4 ...
      + E_pump_defect*Rp_of(Pp_each_CW)*(n_tot-n1_cw-n4_cw);
dT_cw = Qcw*V_xtal/(Hconv*Aside);
fprintf('\n--- CW lumped average temperature, 134 W total pump ---\n');
fprintf('Q_cw = %.3g W/m^3 -> steady dT = %.2f K -> T_avg = %.2f C\n', Qcw, dT_cw, T0_amb-273.15+dT_cw);
fprintf('(paper: ~41-42 C simulated average, 43.8 C measured AT 888NM -- not directly comparable, see header note)\n');

%% ===================== 7. PLOTS 3 & 4: SPATIAL TEMPERATURE MAPS =====================
Nx = 41; Ny = 41;
xs = linspace(-w_x/2, w_x/2, Nx);
ys = linspace(-hgt/2, hgt/2, Ny);
z_near = linspace(0, min(3e-3,lc), 60);        % fine sampling near the entrance face (captures the fast 808nm decay)
z_far  = linspace(z_near(end), lc, 20);        % coarse sampling for the rest of the rod (little happens there)
zs = unique([z_near, z_far(2:end)]);
Nz = numel(zs);
[Xg, Yg] = meshgrid(xs, ys);

Tfull = zeros(Ny, Nx, Nz);

for kz = 1:Nz
    z = zs(kz);
   
    Ip = (2*Pp_each_CW/(pi*omega_p0^2)) * exp(-2*(Xg.^2+Yg.^2)/omega_p0^2) * exp(-alpha0*z);
    Rp_local = sigma_ap0*Ip/(h*nu_p);                          % Eq. 8

    n4_local = solve_n4_local(Rp_local, Ws_cw, n_tot, tau4, tau1, beta4to1, Wup);   % Eq. 2, local steady state
    n1_local = tau1*(beta4to1*n4_local/tau4 + Wup*n4_local.^2 + Ws_cw*n4_local)/(1+Ws_cw*tau1);  % Eq. 5
    nup_local = Wup*n4_local.^2*tauup;                          % Eq. 1
    n0_local = n_tot - n1_local - n4_local;
    Q_local = E10*n1_local/tau1 + E41*nup_local/tauup + E4f*n4_local/tau4 ...
              + E_pump_defect*Rp_local.*n0_local;   % Eq. 18 + 808nm pump-relaxation term (Section 5)

    Tfull(:,:,kz) = solve_2d_conduction(Nx, Ny, w_x, hgt, Kc, Ka, Hconv, T0_amb, Q_local);  % Eq. 19 + Eqs. 21-24
end

fprintf('\n--- Spatial thermal map (134 W total pump / 67 W per crystal) ---\n');
fprintf('Peak temperature = %.2f C (at the pump-entrance face, z=0)\n', max(Tfull(:))-273.15);
fprintf('Volume-averaged temperature = %.2f C\n', mean(Tfull(:))-273.15);
fprintf('(paper: peak ~57-58 C, average ~41-42 C simulated / 43.8 C measured AT 888NM -- not directly comparable, see header note)\n');

save_map(outdir, 'temperature_map_pump_face', xs*1e3, ys*1e3, Tfull(:,:,1)-273.15, ...
         'x [mm]', 'y [mm]', 'Temperature [C]', ...
         'Temperature at pump-entrance face (z=0)');

ix0 = round((Nx+1)/2);
zs_plot = linspace(0, lc, 200);
T_slice = squeeze(Tfull(:,ix0,:));                      % Ny x Nz (graded)
T_slice_plot = interp1(zs, T_slice', zs_plot)';         % Ny x 200 (uniform, for display only)
save_map(outdir, 'temperature_map_longitudinal', zs_plot*1e3, ys*1e3, T_slice_plot-273.15, ...
         'z [mm] (along crystal length, z=0 is pump-entrance face)', 'y [mm]', 'Temperature [C]', ...
         'Longitudinal temperature map at x=0 mid-plane, cf. paper Fig. 3e/f');

%% ===================== 8. PLOTS 5 & 6: THERMO-OPTIC LENSING, Delta n_e MAPS =====================
dne_dT = 8.5e-6;    % [1/K], extraordinary-axis thermo-optic coefficient

n_e_full = ne + dne_dT*(Tfull - T0_amb);   % Koechner Eq. 7.39, applied pointwise to Section 7's T(x,y,z)

dn_e_peak = dne_dT * (max(Tfull(:)) - T0_amb);
fprintf('\n--- Thermo-optic lensing (supplementary; not computed in the paper) ---\n');
fprintf('Peak temperature rise = %.2f K above ambient\n', max(Tfull(:)) - T0_amb);
fprintf('Peak Delta n_e = %.3e   [n_e ranges %.6f to %.6f]\n', dn_e_peak, min(n_e_full(:)), max(n_e_full(:)));

save_map(outdir, 'dn_e_map_pump_face', xs*1e3, ys*1e3, (n_e_full(:,:,1)-ne)*1e6, ...
         'x [mm]', 'y [mm]', '\Delta n_e [\times10^{-6}]', ...
         'Extraordinary-axis \Delta n_e at pump-entrance face (z=0)');

n_e_slice = squeeze(n_e_full(:,ix0,:));                 % Ny x Nz (graded)
n_e_slice_plot = interp1(zs, n_e_slice', zs_plot)';     % Ny x 200 (uniform, for display only)
save_map(outdir, 'dn_e_map_longitudinal', zs_plot*1e3, ys*1e3, (n_e_slice_plot-ne)*1e6, ...
         'z [mm] (along crystal length, z=0 is pump-entrance face)', 'y [mm]', '\Delta n_e [\times10^{-6}]', ...
         '\Delta n_e along the crystal length, x=0 mid-plane');

fprintf('\nSee %s for saved plots/data.\n', outdir);

%% ===================== 9. SMALL-SIGNAL SINGLE-PASS GAIN (no cavity feedback) =====================
% CONTEXT: this answers a specific question from Anya/Chris (gain
% measurement, not the laser cavity above) -- a weak PROBE beam through a
% CW-pumped rod, with no resonator and no circulating cavity photon field
% to saturate the gain. That is the one physical difference from Sections
% 4-8: the stimulated-emission term Ws is set to ZERO (there is no strong
% intracavity field to stimulate emission), leaving a simpler steady-state
% balance -- pump-in vs. spontaneous decay + ETU only:
%   Rp*(ntot-n4) = n4/tau4 + 2*Wup*n4^2        (Eq. 2, with Ws=0)
% which is a quadratic in n4, solved in closed form below. The resulting
% small-signal gain coefficient is:
%   g0 = sigma_es0 * n4        [1/m]   (n1 negligible next to n4, same
%                                        justification as everywhere else
%                                        in this script, Section 4)
% DOPING: Anya's crystal is 1 at.% Nd:YVO4, not the paper's 0.5 at.% --
% ion density n_tot scales linearly with doping fraction (it's a
% substitutional dopant on a fixed lattice site density), so
% ntot_gain = n_tot * (doping_pct/0.5). sigma_ap0, sigma_es0, tau4, Wup
% are intrinsic per-ion properties and do NOT scale with doping.
% NOTE: this section automatically reflects the 808nm sigma_ap0/lambda_p
% set in Sections 2/6 above -- no separate wavelength choice is made here.
doping_pct   = 1.0;      % at.% doping for this measurement
Pabs_gain    = 10;       % W, absorbed pump power
omega_p0_gg  = 450e-6;   % m, pump waist -- CONFIRM against Anya's actual setup; g0 scales ~1/w^2

ntot_gain = n_tot * (doping_pct/0.5);
I0_gain   = 2*Pabs_gain/(pi*omega_p0_gg^2);      % peak on-axis intensity, same Gaussian form as Section 7
Rp_gain   = sigma_ap0*I0_gain/(h*nu_p);          % Eq. 8

a_q = 2*Wup;
b_q = 1/tau4 + Rp_gain;
c_q = -Rp_gain*ntot_gain;
n4_gain = (-b_q + sqrt(b_q^2 - 4*a_q*c_q))/(2*a_q);   % Eq. 2, Ws=0, solved for n4
g0 = sigma_es0*n4_gain;                                % [1/m]

fprintf('\n--- Small-signal single-pass gain (no cavity feedback) ---\n');
fprintf('Doping = %.2f%%, n_tot = %.3e m^-3, Pabs = %.1f W, w_p0 = %.0f um\n', ...
        doping_pct, ntot_gain, Pabs_gain, omega_p0_gg*1e6);
fprintf('I0 = %.1f W/cm^2, Rp = %.2f 1/s\n', I0_gain/1e4, Rp_gain);
fprintf('n4 = %.3e m^-3 (%.2f%% of ions excited)\n', n4_gain, n4_gain/ntot_gain*100);
fprintf('g0 = %.1f 1/m  =  %.2f %%/mm  =  %.1fx over 1 cm\n', ...
        g0, g0*1e-3*100, exp(g0*1e-2));

%% ===================== 10. THRESHOLD PUMP INTENSITY FOR A TARGET GAIN (inverse of Section 9) =====================
% Same physics as Section 9, run backwards: given a target small-signal
% gain coefficient (as %/mm), solve for the pump INTENSITY that produces
% it, at the same doping. Inverting Rp*(ntot-n4) = n4/tau4 + 2*Wup*n4^2
% (Eq. 2, Ws=0) for Rp given a target n4 = g0_target/sigma_es0:
target_pct_per_mm = 5.0;    % the gain threshold in question, e.g. "gain > 5%" per mm

g0_target = target_pct_per_mm * 10;         % %/mm -> 1/m  (g0[1/m]*1e-3*100 = %/mm)
n4_target = g0_target/sigma_es0;
Rp_target = (n4_target/tau4 + 2*Wup*n4_target^2) / (ntot_gain - n4_target);
I0_target = Rp_target*h*nu_p/sigma_ap0;      % Eq. 8, inverted for intensity

fprintf('\n--- Threshold pump intensity for g0 > %.1f %%/mm (%.1f%% doping) ---\n', ...
        target_pct_per_mm, doping_pct);
fprintf('Required n4 = %.3e m^-3 (%.4f%% of ions), Rp = %.2f 1/s\n', ...
        n4_target, n4_target/ntot_gain*100, Rp_target);
fprintf('Threshold pump intensity I0 = %.3e W/m^2  =  %.1f W/cm^2\n', I0_target, I0_target/1e4);
fprintf('(equivalent absorbed power at w_p0 = %.0f um waist: %.2f W)\n', ...
        omega_p0_gg*1e6, I0_target*pi*omega_p0_gg^2/2);

%% ===================== Plot 7: gain coefficient vs. pump intensity =====================
% Direct visual answer to "what pump intensity do we start to see gain
% > X%?" -- sweeps I0 across a range spanning well below/above I0_target
% (Section 10) and solves the same Eq. 2 (Ws=0) balance at each point.
I0_sweep = linspace(0.05*I0_target, 3*I0_target, 200);
g0_sweep = zeros(size(I0_sweep));
for k = 1:numel(I0_sweep)
    Rp_k = sigma_ap0*I0_sweep(k)/(h*nu_p);
    b_k = 1/tau4 + Rp_k;
    c_k = -Rp_k*ntot_gain;
    n4_k = (-b_k + sqrt(b_k^2 - 4*a_q*c_k))/(2*a_q);
    g0_sweep(k) = sigma_es0*n4_k*1e-3*100;   % [1/m] -> %/mm
end

f7 = figure('Name', 'Small-signal gain vs. pump intensity');
plot(I0_sweep/1e4, g0_sweep, 'o-', 'LineWidth', 1.5); hold on; grid on;
yline(target_pct_per_mm, 'r--', sprintf('%.1f %%/mm target', target_pct_per_mm), 'LineWidth', 1.2);
xline(I0_target/1e4, 'r--', sprintf('%.0f W/cm^2', I0_target/1e4), 'LineWidth', 1.2);
plot(I0_target/1e4, target_pct_per_mm, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
xlabel('Pump intensity I_0 [W/cm^2]'); ylabel('Small-signal gain coefficient g_0 [%/mm]');
title(sprintf('Gain vs. pump intensity (%.1f%% doping) -- crosses %.1f%%/mm at %.0f W/cm^2', ...
      doping_pct, target_pct_per_mm, I0_target/1e4));
pngfile7 = fullfile(outdir, 'gain_vs_intensity.png');
print(f7, pngfile7, '-dpng', '-r150');
fprintf('Displayed and saved plot: %s\n', pngfile7);

%% ===================== LOCAL FUNCTIONS =====================

function n1 = n1_qss(n4, Ws, p)
% Adiabatic quasi-steady-state value of n1 (Eq. 5). 
    n1 = p.tau1 .* (p.beta4to1.*n4./p.tau4 + p.Wup.*n4.^2 + Ws.*n4) ./ (1 + Ws.*p.tau1);
end

function [Ps, phi, n4, n1] = cw_state(Pp_each, Rp_of, Dn_th, n_tot, tau4, tau1, ...
                                       beta4to1, Wup, nc, c, sigma_es0, ...
                                       nu_s, Vm, tau_r, Toc)
% Closed-form CW steady state
    h = 6.62607015e-34;
    Rp = Rp_of(Pp_each);
    n0_approx = n_tot - Dn_th;
    numer = Rp*n0_approx - Dn_th/tau4 - 2*Wup*Dn_th^2;
    if numer <= 0
        % below threshold: negligible stimulated emission in this
        % simplified CW treatment
        Ps = 0; phi = 0; n4 = Rp*n0_approx*tau4; n1 = 0;
        return;
    end
    Ws = numer/Dn_th;                              % Eq. 2, solved for Ws
    phi = Ws/((c/nc)*sigma_es0);                    % Eq. 9
    n4 = Dn_th;
    n1 = tau1*(beta4to1*n4/tau4 + Wup*n4^2 + Ws*n4)/(1+Ws*tau1);   % Eq. 5
    Ps = h*nu_s*phi*Vm/tau_r*log(1/(1-Toc));        % Eq. 16
end

function n4 = solve_n4_local(Rp, Ws, n_tot, tau4, tau1, beta4to1, Wup)
    n1fun = @(n4) tau1*(beta4to1*n4/tau4 + Wup*n4.^2 + Ws*n4)./(1+Ws*tau1);   % Eq. 5

    n4 = Rp*n_tot*tau4 ./ (1 + Rp*tau4 + Ws*tau4);   % initial guess (linearized, ignoring ETU/n1)
    for iter = 1:25
        n1 = n1fun(n4);
        eps_ = n4*1e-6 + 1e-10;
        dn1 = (n1fun(n4+eps_) - n1) ./ eps_;          % numerical derivative, vectorized
        n0 = n_tot - n1 - n4;
        f  = Rp.*n0 - n4/tau4 - 2*Wup*n4.^2 - Ws*(n4-n1);   % Eq. 2
        dfdn4 = Rp.*(-dn1-1) - 1/tau4 - 4*Wup*n4 - Ws*(1-dn1);
        n4 = n4 - f./dfdn4;
        n4 = min(max(n4, 0), n_tot);
    end
end

function T = solve_2d_conduction(Nx, Ny, Lx, Ly, Kx, Ky, H, T0, Qgrid)
    dx = Lx/(Nx-1); dy = Ly/(Ny-1);
    N = Nx*Ny;
    idx = @(i,j) (j-1)*Nx + i;    % 1-based linear index, i: x, j: y

    maxnnz = 5*N;
    rows = zeros(maxnnz,1); cols = zeros(maxnnz,1); vals = zeros(maxnnz,1);
    b = zeros(N,1);
    p = 0;

    for j = 1:Ny
        for i = 1:Nx
            k = idx(i,j);
            diagv = 0;
            if i>1
                g = Kx*dy/dx; p=p+1; rows(p)=k; cols(p)=idx(i-1,j); vals(p)=-g; diagv=diagv+g;
            else
                g = H*dy; diagv=diagv+g; b(k)=b(k)+g*T0;
            end
            if i<Nx
                g = Kx*dy/dx; p=p+1; rows(p)=k; cols(p)=idx(i+1,j); vals(p)=-g; diagv=diagv+g;
            else
                g = H*dy; diagv=diagv+g; b(k)=b(k)+g*T0;
            end
            if j>1
                g = Ky*dx/dy; p=p+1; rows(p)=k; cols(p)=idx(i,j-1); vals(p)=-g; diagv=diagv+g;
            else
                g = H*dx; diagv=diagv+g; b(k)=b(k)+g*T0;
            end
            if j<Ny
                g = Ky*dx/dy; p=p+1; rows(p)=k; cols(p)=idx(i,j+1); vals(p)=-g; diagv=diagv+g;
            else
                g = H*dx; diagv=diagv+g; b(k)=b(k)+g*T0;
            end
            p=p+1; rows(p)=k; cols(p)=k; vals(p)=diagv;
            b(k) = b(k) + Qgrid(j,i)*dx*dy;
        end
    end
    rows = rows(1:p); cols = cols(1:p); vals = vals(1:p);
    A = sparse(rows, cols, vals, N, N);
    Tvec = A\b;
    T = reshape(Tvec, Nx, Ny)';   % back to (Ny x Nx), matching Qgrid's orientation
end

function save_plot(outdir, name, x, y, xlab, ylab, titlestr)
% Renders a normal visible figure window AND saves a PNG copy to disk;
% falls back to CSV only if this MATLAB has no raster backend available
% (e.g. some restricted MATLAB Online configurations).
    pngfile = fullfile(outdir, [name '.png']);
    csvfile = fullfile(outdir, [name '.csv']);
    saved_png = false;
    try
        f = figure('Name', titlestr);
        plot(x, y, 'o-', 'LineWidth', 1.5); grid on;
        xlabel(xlab); ylabel(ylab); title(titlestr);
        print(f, pngfile, '-dpng', '-r150');
        saved_png = true;
        fprintf('Displayed and saved plot: %s\n', pngfile);
    catch ME
        fprintf('Could not render "%s" as a PNG in this environment (%s).\n', name, ME.message);
    end
    if ~saved_png
        T = table(x(:), y(:), 'VariableNames', {matlab.lang.makeValidName(xlab), matlab.lang.makeValidName(ylab)});
        try
            writetable(T, csvfile);
            fprintf('Saved data instead as CSV: %s  (%s vs %s)\n', csvfile, xlab, ylab);
        catch
            fprintf('%s , %s\n', xlab, ylab);
            for i = 1:numel(x)
                fprintf('%.6g , %.6g\n', x(i), y(i));
            end
        end
    end
end

function save_map(outdir, name, xvec, yvec, Zdata, xlab, ylab, cbarlab, titlestr)
% Same as save_plot, for a 2D heatmap (imagesc-style).
    pngfile = fullfile(outdir, [name '.png']);
    csvfile = fullfile(outdir, [name '.csv']);
    saved_png = false;
    try
        f = figure('Name', titlestr);
        imagesc(xvec, yvec, Zdata); axis image; set(gca,'YDir','normal');
        cb = colorbar; cb.Label.String = cbarlab;
        xlabel(xlab); ylabel(ylab); title(titlestr);
        print(f, pngfile, '-dpng', '-r150');
        saved_png = true;
        fprintf('Displayed and saved map: %s\n', pngfile);
    catch ME
        fprintf('Could not render "%s" as a PNG in this environment (%s).\n', name, ME.message);
    end
    if ~saved_png
        M = [0, xvec(:)'; yvec(:), Zdata];
        try
            writematrix(M, csvfile);
            fprintf('Saved data instead as CSV: %s  (row1=x, col1=y, corner=0)\n', csvfile);
        catch
            fprintf('Could not save "%s" as CSV either; skipping.\n', name);
        end
    end
end

function E = wn2J(wn_cm1, h, c)
% Convert a wavenumber in cm^-1 to a photon energy in Joules.
    E = h*c*(wn_cm1*100);
end