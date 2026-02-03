function Analyze_Stability(L, Rc1, Rc2)
    % Stability
    g1 = 1 - L/Rc1; % stability parameter 1
    g2 = 1 - L/Rc2; % stability paramter 2
    g = g1*g2; % stability product, 0 < g < 1 for a stable cavity

    if g > 0 && g < 1
        disp('The cavity is stable, g = %.2f' , g);
    else
        disp('Unstable cavity! g = %.2f', g)
    end
end