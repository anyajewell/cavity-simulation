function [gain_medium] = Initialize_Gain_Medium(sim, mirror)
    g0 = 3;
    I_sat = 10^8; % gain medium saturation intensity, [W/m^2]
    sigma = mirror(1).D/4;
    g0_profile = g0*exp(-(sim.X.^2+sim.Y.^2)./(2*sigma^2)); % gain function
    Aeff = 2*pi*sigma^2; % effective area of gain medium
    P_sat = I_sat*Aeff; % saturation power
    gain_medium.I_sat = I_sat; gain_medium.g0_profile = g0_profile; gain_medium.P_sat = P_sat;
end
    

