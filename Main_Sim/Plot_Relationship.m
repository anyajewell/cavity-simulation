function fig = Plot_Relationship(var1, var2)
    fig = figure;
    plot(var1, var2, 'LineWidth', 1.5)
    title({'Fresnel number versus' 'converged loss'});
    xlabel('N_F');
    ylabel('Loss per round trip');
    grid on;
end