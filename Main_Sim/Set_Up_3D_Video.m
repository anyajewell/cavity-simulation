function [video, frameStride, fig, layout, sliceSurfaces, ix, iy, sliceIndices, intensityFactor] = Set_Up_3D_Video(sim, consts, laser)

    % Transverse crop around the optical axis
    ix = find(abs(sim.X(1,:)) <= 1);
    iy = find(abs(sim.Y(:,1)) <= 1);
    x = sim.X(1,ix);
    y = sim.Y(iy,1);
    [Xslice, Yslice] = meshgrid(x, y);
    
    % Set up however many slices you want to see
    z = sim.z(:);
    zTargets = linspace(z(1), z(end), 10);
    [~, sliceIndices] = min(abs(z - zTargets), [], 1);
    sliceIndices = unique(sliceIndices, 'stable');
    zSlice_km = z(sliceIndices) / 1e3;
    zLimits_km = ([z(1), z(end)] + [-0.5, 0.5]*sim.dz) / 1e3;
    
    intensityFactor = 0.5 * consts.c * consts.eps0;
    Iplus = intensityFactor * abs(laser.Eplus(iy,ix,sliceIndices)).^2;
    Iminus = intensityFactor * abs(laser.Eminus(iy,ix,sliceIndices)).^2;
    
    % Use the same fixed color scale for both directions and all frames.
    intensityMax = max(max(Iplus, [], 'all'), max(Iminus, [], 'all'));
    if intensityMax == 0
        intensityMax = 1;
    end
    
    frameStride = 1;
    fig = figure('Color', 'w', 'Position', [100 100 1400 750]);
    layout = tiledlayout(fig, 2, 1, ...
        'TileSpacing', 'compact', 'Padding', 'compact');
    colormap(fig, parula);
    
    sliceSurfaces = gobjects(2, numel(sliceIndices));
    directionLabels = {'Right-going field', 'Left-going field'};
    intensity = {Iplus, Iminus};
    
    for d = 1:2
        ax = nexttile(layout);
        hold(ax, 'on');
        volume = intensity{d};
    
        for j = 1:numel(sliceIndices)
            % Plot coordinates are (cavity z, transverse x, transverse y).
            Zplane = zSlice_km(j) * ones(size(Xslice));
            sliceSurfaces(d,j) = surf(ax, Zplane, Xslice, Yslice, ...
                volume(:,:,j), 'FaceColor', 'texturemap', ...
                'EdgeColor', 'none');
        end
    
        xlabel(ax, 'Cavity position z [km]');
        ylabel(ax, 'Transverse x [m]');
        zlabel(ax, 'Transverse y [m]');
        title(ax, directionLabels{d});
    
        xlim(ax, zLimits_km);
        ylim(ax, [min(x), max(x)]);
        zlim(ax, [min(y), max(y)]);
        xticks(ax, zSlice_km);
        clim(ax, [0, intensityMax]);
    
        pbaspect(ax, [4 1 1]);
        view(ax, [-35 22]);
        grid(ax, 'on');
        box(ax, 'on');

        % Black tick labels, axis labels, and panel title
        set(ax, 'XColor', 'k', 'YColor', 'k', 'ZColor', 'k');
        ax.XLabel.Color = 'k';
        ax.YLabel.Color = 'k';
        ax.ZLabel.Color = 'k';
        ax.Title.Color = 'k';

        % White grid lines against the dark cavity walls
        ax.Color = 'w'; % White floor and walls
        ax.GridColor = [0.82 0.82 0.82]; % Light gray lines visible on white
        ax.GridAlpha = 1;
    end
    
    cb = colorbar(ax);
    cb.Layout.Tile = 'east';
    cb.Label.String = 'Intensity [W/m^2]';
    cb.Color = 'k';
    cb.Label.Color = 'k';
    
    video = VideoWriter('cavity_3D_evolution.mp4', 'MPEG-4');
    video.FrameRate = 24;
    open(video);

end