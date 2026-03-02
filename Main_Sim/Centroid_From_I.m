function [cx, cy] = Centroid_From_I(sim, I)
    % I is size [Ny, Nx] same as imagesc(sim.x, sim.y, I)
    
    P  = sum(I(:));
    Ix = sum(I, 1);          % 1 x Nx   (sum over y)
    Iy = sum(I, 2);          % Ny x 1   (sum over x)
    
    cx = sum(sim.x(:)'.*Ix) / P;
    cy = sum(sim.y(:) .*Iy) / P;
end