function fig = Plot_Center(center, zs)
    % Plotting
    fig = figure;
    plot(center,zs,'b','LineWidth',1.5); hold on;
    xlabel('Displacement [m]'); ylabel('Distance [m]');
    title('Beam Center Over Propagation')
end