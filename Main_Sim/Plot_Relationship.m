function fig = Plot_Relationship(var1, var2)
    fig = figure;
    plot(var1, var2, '.', 'MarkerSize', 15)
    title({'Gain medium profile width' 'versus converged loss fraction'});
    xlabel('\sigma / R');
    ylabel('Converged loss fraction');
    grid on;
end