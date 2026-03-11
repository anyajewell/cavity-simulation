function fig = Plot_Relationship(var1, var2)
    fig = figure;
    plot(var1, var2, '.', 'MarkerSize', 15); hold on
    %plot(var1, var3, '.', 'MarkerSize', 15, 'DisplayName', 'mirror 2 loss');
    title({'Gain medium profile width' 'versus converged round-trip loss fraction'});
    xlabel('\sigma / R');
    ylabel('Converged loss fraction');
    %legend('Location','best')
    grid on;
end