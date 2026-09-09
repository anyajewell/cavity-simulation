function fig = Plot_Relationship(var1, var2)
    fig = figure;
    plot(var1, var2, 'LineWidth', 1.5)
    title({'Fresnel Number vs.' 'Converged Loss'});
    xlabel('N_F');
    ylabel('Loss per round trip');
    grid on;
end