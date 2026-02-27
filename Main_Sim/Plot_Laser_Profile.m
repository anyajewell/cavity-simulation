function fig = Plot_Laser_Profile(sim, laser)
    
    fig = figure;
    set(fig, 'Color', 'w');
    set(gca, 'Color', 'w');
    imagesc(sim.x, sim.y, abs(laser.Gau));
    axis equal
    axis tight
    axis([-0.5 0.5 -0.5 0.5]);
    colorbar
    xlabel('x [m]')
    ylabel('y [m]')
    title('|E(x,y)|')

end