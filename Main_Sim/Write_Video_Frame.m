function outputs = Write_Video_Frame(sim, laser, toggles, outputs)
    % Visualization
    imagesc(sim.x,sim.y,abs(laser.Gau)); axis([-0.5 0.5 -0.5 0.5]); axis square; xlabel('x [m]'); ylabel('y [m]'); hold on; 
    if toggles.track_centers == true
        plot(outputs.centerx(end),outputs.centery(end),'ro'); hold off;
    end
    f = getframe(gcf); display(laser.pos);
    writeVideo(outputs.v,f);
end