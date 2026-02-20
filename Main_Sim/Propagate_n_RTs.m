function [laser, outputs, gain_medium] = Propagate_n_RTs(consts, sim, laser, frame, mirror, outputs, toggles, gain_medium)

    for a = 1:sim.RTs
        Gau_a = laser.Gau; % beam profile at the start of this round trip

        if sign(sim.dz) > 0
            [laser, outputs] = R(consts, sim, laser, frame, mirror, outputs, toggles, gain_medium, a);
            % Visualization
            imagesc(sim.x,sim.y,abs(laser.Gau)); axis([-0.5 0.5 -0.5 0.5]); axis square; xlabel('x [m]'); ylabel('y [m]'); hold on; 
            if toggles.track_centers == true
                plot(outputs.centerx(end),outputs.centery(end),'ro'); hold off;
            end
            f = getframe(gcf); display(sim.z(i));
            writeVideo(outputs.v,f);
            
        else
            L()

        end
        

        % Visualization
        imagesc(sim.x,sim.y,abs(laser.Gau)); axis([-0.5 0.5 -0.5 0.5]); axis square; xlabel('x [m]'); ylabel('y [m]');  
        set(gcf, 'Color', 'w'); set(gca, 'Color', 'w'); hold on;
        if toggles.track_centers == true
            plot(outputs.centerx(end),outputs.centery(end),'ro'); hold off;
        end
        f = getframe(gcf); display(sim.z(i));
        writeVideo(outputs.v,f);
   
        % Calculate and store loss
        loss_a =  1 - sum(abs(laser.Gau).^2,'all') / sum(abs(Gau_a).^2,'all'); % loss fraction this round-trip
        outputs.loss_frac(a) = loss_a;

    end

    if toggles.track_centers == true % fix counting
        outputs.centerx(end) = [];
        outputs.centery(end) = [];
    end

end