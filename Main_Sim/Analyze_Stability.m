function Analyze_Stability(sim, mirror)
    % Stability
    g1 = 1 - sim.L/mirror(1).Rc; % stability parameter 1
    g2 = 1 - sim.L/mirror(2).Rc; % stability paramter 2
    g = g1*g2; % stability product, 0 < g < 1 for a stable cavity

    if g > 0 && g < 1
        fprintf('The cavity is stable, g = %.2f', g);
    elseif g == 0 || g == 1
        fprintf('Marginally stable, g = %.2f', g)
    else
        fprintf({'Unstable cavity! g = %.2f', g})
    end
end