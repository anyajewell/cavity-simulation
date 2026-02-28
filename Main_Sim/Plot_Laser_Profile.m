function fig = Plot_Laser_Profile(sim, profile)
    
    fig = figure;
    set(fig, 'Color', 'w');
    set(gca, 'Color', 'w');
    imagesc(sim.x, sim.y, abs(profile).^2);
    axis equal
    axis tight
    axis([-0.5 0.5 -0.5 0.5]);
    colorbar
    xlabel('x [m]')
    ylabel('y [m]')
    title('I(x,y)')

end