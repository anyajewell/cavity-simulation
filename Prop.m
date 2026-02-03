function [Gau, centerx, centery] = Prop(Gau, Omega, accel, dt, c, Ld, dx, dz, x, y, z, X, Y, k0, v, centerx, centery, i)
    
    %v_perp = v_perp + accel*dt; % transverse velocity
    %dx_frame = v_perp*dt; % frame translation this step
    dv = accel*dt; % step in velocity
    dTheta = dv/c; % tilt due to transverse acceleration

    % Propagation using angular spectrum method (Schmidt et al.)
    [x2dum, Gau] = ang_spec_prop(Gau,Ld,dx,dx,dz);
    
    % New X-position after trajectory along characteristic curve
    rot_shift = 0.5*Omega/c * (z(i+1).^2-z(i).^2);
    Xnew = X - rot_shift;
    
    % Analytic integration of phase term along characteristic 
    Xint = X*(z(i+1)-z(i)) - 0.5*Omega/c * ...
        (1/3 * (z(i+1).^3-z(i).^3) + 0.5*z(i)*(z(i+1)-z(i)));
    
    % Phase screen implementation for tilt terms
    A = exp(-1i * k0 * Omega / c * Xint); % rotational tilt phase screen
    B = exp(-1i * k0 * dTheta.*X); % transverse acceleration tilt phase screen
    Gau = Gau .* A .* B; % apply tilt
    
    % Implementation of shift interpolation
    Gau = interp2(Xnew,Y,Gau,X,Y,'spline'); % interpolate onto the new grid
    
    centerx(i+1) = trapz(trapz(X.*abs(Gau).^2))/trapz(trapz(abs(Gau).^2)); % track center x
    centery(i+1) = trapz(trapz(Y.*abs(Gau).^2))/trapz(trapz(abs(Gau).^2)); % track center y

    imagesc(x,y,abs(Gau)); axis([-0.5 0.5 -0.5 0.5]); axis square; xlabel('x [m]'); ylabel('y [m]'); 
    hold on; plot(centerx(i+1),centery(i+1),'ro'); hold off; frame = getframe(gcf); display(z(i));
    writeVideo(v,frame);
      
end

