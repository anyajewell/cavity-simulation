function fig = Plot_Results_vs_Analytic(centerx, x_ana, z)
    % Plotting
    fig = figure;
    title({'Displacement in x Over Propagation'})
    plot(centerx,z,'b','DisplayName','Numeric'); hold on;
    plot(x_ana,z,'r','DisplayName','Analytic')
    xlabel('Displacement [m]'); ylabel('Distance [m]');
    legend();
end