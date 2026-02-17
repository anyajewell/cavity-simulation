function [laser, outputs] = Prop(consts, sim, laser, frame, outputs, toggles, i)

    dv = frame.accel*sim.dt; % step in velocity
    dTheta = dv/consts.c; % tilt due to transverse acceleration

    % Propagation using angular spectrum method (Schmidt et al.)
    [~, laser.Gau] = ang_spec_prop(laser.Gau,laser.Ld,sim.dx,sim.dx,sim.dz);
    
    % New X-position after trajectory along characteristic curve
    rot_shift = 0.5*frame.Omega/consts.c * (sim.z(i+1).^2-sim.z(i).^2);
    Xnew = sim.X - rot_shift;
    %xnew = sim.x - rot_shift;
    
    % Analytic integration of phase term along characteristic 
    Xint = sim.X*(sim.z(i+1)-sim.z(i)) - 0.5*frame.Omega/consts.c * ...
        (1/3 * (sim.z(i+1).^3-sim.z(i).^3) + 0.5*sim.z(i)*(sim.z(i+1)-sim.z(i)));
    
    % Phase screen implementation for tilt terms
    A = exp(-1i * laser.k0 * frame.Omega / consts.c * Xint); % rotational tilt phase screen
    B = exp(-1i * laser.k0 * dTheta.*sim.X); % transverse acceleration tilt phase screen
    laser.Gau = laser.Gau .* A .* B; % apply tilt
    
    % Implementation of shift interpolation
    %laser.Gau = interp1(xnew, laser.Gau.', sim.x, 'linear', 0).';
    laser.Gau = interp2(Xnew,sim.Y,laser.Gau,sim.X,sim.Y,'spline'); % interpolate onto the new grid
    
    if toggles.track_centers == true
        outputs.centerx(end+1) = trapz(trapz(sim.X.*abs(laser.Gau).^2))/trapz(trapz(abs(laser.Gau).^2)); % track center x
        outputs.centery(end+1) = trapz(trapz(sim.Y.*abs(laser.Gau).^2))/trapz(trapz(abs(laser.Gau).^2)); % track center y
    end

    % Plot every step here
    % imagesc(x,y,abs(Gau)); axis([-0.5 0.5 -0.5 0.5]); axis square; xlabel('x [m]'); ylabel('y [m]'); 
    % set(gcf, 'Color', 'w'); set(gca, 'Color', 'w'); hold on; 
    % if track_centers == true
    %     plot(centerx(end),centery(end),'ro'); hold off;
    % end
    % 
    % frame = getframe(gcf); display(z(i));
    % writeVideo(v,frame);

end

