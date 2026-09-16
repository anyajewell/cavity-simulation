function E = Apply_Thermo_Optic_Effect(E)

    OPD_th = sum(delta_n,3) * dz;
    
    % Remove piston phase
    OPD_th = OPD_th - OPD_th(center_idx,center_idx);
    
    phi_th = k0 * OPD_th;
    
    E = E .* exp(1i*phi_th);

end
