function Write_Cavity_Frame(video, fig, layout, sliceSurfaces, ...
    laser, ix, iy, sliceIndices, intensityFactor, t)

    Iplus  = intensityFactor * abs(laser.Eplus(iy,ix,sliceIndices)).^2;
    Iminus = intensityFactor * abs(laser.Eminus(iy,ix,sliceIndices)).^2;
    
    for j = 1:numel(sliceIndices)
        sliceSurfaces(1,j).CData = Iplus(:,:,j);
        sliceSurfaces(2,j).CData = Iminus(:,:,j);
    end
    
    title(layout, sprintf('Transverse intensity | t = %.3f ms', 1e3*t), ...
        'Color', 'k');
    drawnow;
    writeVideo(video, getframe(fig));

end