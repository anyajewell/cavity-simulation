function fig = Plot_Relationship(var1, var2)

    % Linear fit to measured data
    p = polyfit(var1, var2, 1);

    % Extrapolate to 5.5 A
    I_ext = linspace(max(var1), 5.5, 100);
    P_ext = polyval(p, I_ext);

    % Plot
    fig = figure;
    hold on;

    % Measured data
    plot(var1, var2, 'b-', 'LineWidth', 1.5);

    % Extrapolated data
    plot(I_ext, P_ext, 'r--', 'LineWidth', 1.5);

    title('Current vs. Power for 808 Diode');
    xlabel('I [A]');
    ylabel('P [W]');

    legend('Measured', 'Linear Extrapolation', 'Location', 'best');

    xlim([min(var1), 5.5]);
    grid on;

    % Report extrapolated power
    fprintf('Predicted power at 5.5 A: %.2f W\n', polyval(p, 5.5));

end
