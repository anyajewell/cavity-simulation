function fig = Plot_Relationship(var1, var2, var3)
    fig = figure;
    plot(var1, var2, '.', 'MarkerSize', 15, 'DisplayName', 'mirror 1 loss'); hold on
    plot(var1, var3, '.', 'MarkerSize', 15, 'DisplayName', 'mirror 2 loss');
    title({'Gain medium profile width' 'versus converged loss fraction'});
    xlabel('\sigma / R');
    ylabel('Converged loss fraction');
    legend('Location','best')
    grid on;
end