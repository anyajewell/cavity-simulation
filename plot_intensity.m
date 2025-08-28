function fig = plot_intensity(E, consts, Z_traveled, Z_position, step, X, Y, x_circ2, y_circ2)

    fig = figure; 
    I = 0.5*consts.eps0*consts.c*abs(E).^2;
    surf(X, Y, I, 'LineStyle','none');
    hold on;
    plot3(x_circ2, y_circ2, zeros(size(x_circ2)), 'r-', 'LineWidth', 2); % plot mirror outline
    hold off;
    title({
        sprintf('Laser Mode at Z = %.1f m', Z_traveled(step)), ...
        sprintf('Intra-Cavity Position = %.1f', Z_position(step))
    })
    zlabel('Intensity')
    xlabel('X [m]')
    ylabel('Y [m]')

end