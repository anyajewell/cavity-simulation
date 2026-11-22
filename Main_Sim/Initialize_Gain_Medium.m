function [gain_medium] = Initialize_Gain_Medium(sim, mirror, gain_type)
    
    if strcmp(gain_type, 'basic') % basic gain profile

        g0 = 3;
        I_sat = 10^8; % gain medium saturation intensity, [W/m^2]
        R = mirror(1).D/2;
        sigma_R = 0.5; % sigma / R
        sigma = R*sigma_R;
        g0_profile = g0*exp(-(sim.X.^2+sim.Y.^2)./(2*sigma^2)); % gain function
        Aeff = 2*pi*sigma^2; % effective area of gain medium
        P_sat = I_sat*Aeff; % saturation power
        gain_medium.I_sat = I_sat; gain_medium.g0_profile = g0_profile; gain_medium.P_sat = P_sat; gain_medium.sigma = sigma;

    else % sophisticated gain model

        gain_medium.R = 1.5e-3; % crystal radius, [m]                     
        gain_medium.L_c = 15e-3; % crystal length, [m] 
        gain_medium.N = 1.26e26; % Nd density at 1 at%, [1/m^3]
        gain_medium.kc = 12.1; % k along c (x) at 25 C, W/(m K) [SOURCED] [T] Table 2
        gain_medium.ka = 8.9; % k along a (y, z) at 25 C, W/(m K) [SOURCED] [T] Table 2
        gain_medium.T_ref = 298.15; % temperature of kc, ka, K [SOURCED] [T] (25 C)
        gain_medium.n_k = 1; % k(T) exponent, scalar or [n_c n_a n_a] [ESTIMATE] phonon-limited k ~ 1/T
        gain_medium.dne_dT = 2.9e-6; % dn_e/dT, [1/K] [SOURCED] [D]; peer-reviewed 633 nm values are 9.5-12.5e-6

        % Copper mount: trapezoidal block, chiller on the top face, crystal centered in a thin contact layer
        gain_medium.W_top = 12.7e-3; % top width, [m] [PLACEHOLDER]
        gain_medium.W_bot = 6.35e-3; % bottom width, [m] [PLACEHOLDER]
        gain_medium.H_blk = 12.7e-3; % block height, [m] [PLACEHOLDER]
        gain_medium.k_cu = 400; % copper conductivity, W/(m K) [SOURCED] standard value
        gain_medium.h_c = 2e4; % contact conductance, W/(m^2 K) [PLACEHOLDER] indium foil ~1e4-5e4
        gain_medium.t_c = 0.1e-3; % contact layer thickness, m (numerical device; see boundary conditions above)
        gain_medium.T0 = 16; % chiller temperature, C                 

        % Pump laser
        gain_medium.P = 10; % pump power at the crystal, [W] 
        gain_medium.lambda_p = 808e-9; % pump wavelength, [m] 
        gain_medium.alpha = 2300; % low-intensity absorption, [1/m] [SOURCED] [D] pi-polarized peak (upper bound)
        gain_medium.w0 = 0.75e-3; % pump waist radius, [m] 
        gain_medium.z_w = L_c/2; % waist position from the pump face, [m]
        gain_medium.M2_p = 44.9; % pump beam quality                
        gain_medium.n_p = 2.1858; % n_e at 808 nm [SOURCED] LaserComponents Nd:YVO4 datasheet

        % Spectroscopy
        gain_medium.tau_f = 90e-6; % fluorescence lifetime, [s] [SOURCED] [D]
        gain_medium.tau_rad = 104.29e-6; % radiative lifetime, [s] [ESTIMATE] low-doping lifetime (Shen Table 1),
        gain_medium.lambda_f = 1032e-9; % mean fluorescence wavelength, [m] [ESTIMATE] Shen Table 1
        gain_medium.W_up = 5.0e-22; % ETU coefficient, [m^3/s] [SOURCED] [C] 0.97 at%, room temperature

        % Laser (lasing = false gives the non-lasing case)
        lasing = true;
        lambda_L = 1064.3e-9; % laser wavelength, [m] [SOURCED] [D]
        sigma_L = 15.6e-23;  % emission cross section (pi), [m^2] [SOURCED] Sato & Taira 2012, as quoted in Shen Table 1
        T_oc = 0.10; % output coupler transmission [CHOSEN]
        L_i = 0.02; % round-trip internal loss [ESTIMATE]
        w_L = 0.8*w0; % laser mode radius in the crystal, [m] [ESTIMATE] design ratio 0.7-1.0 of the pump radius

    end

end
    

