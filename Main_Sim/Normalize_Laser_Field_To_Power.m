function [E_scaled] = Normalize_Laser_Field_To_Power(E, Ptarget, dx, dy, c, eps0)
    I = 0.5*c*eps0*abs(E).^2; % [W/m^2] if E is [V/m]
    Pcur = dx*dy*sum(I,'all'); % current power, [W]
    E_scaled = E*sqrt(Ptarget / Pcur); % rescale E field
end