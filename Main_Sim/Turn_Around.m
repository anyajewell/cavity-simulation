function [Zmax, dz, z] = Turn_Around(Zmax, dz, Nz)

    % Turn around
    Zmax = -Zmax;
    Z0 = -Zmax;
    dz = -dz;
    z = linspace(Z0,Zmax,Nz);

end