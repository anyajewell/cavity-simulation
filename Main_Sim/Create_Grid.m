function [x, y, X, Y] = Create_Grid(N, Lx, dx)    
    % Create grid
    for n=1:N
        for m=1:N
            x(m)=(m-1)*dx-Lx/2;
            y(n)=(n-1)*dx-Lx/2;
        end
    end
    [X,Y] = meshgrid(x,y);
end