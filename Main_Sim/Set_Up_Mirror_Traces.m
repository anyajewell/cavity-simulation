function [x_circ1, y_circ1, x_circ2, y_circ2] = Set_Up_Mirror_Traces(D1, D2)
    % Set up mirror circles for plotting
    r1 = D1/2; % radius of mirror 1
    r2 = D2/2; % radius of mirror 2
    theta = linspace(0,2*pi,400);
    x_circ1 = r1*cos(theta); y_circ1 = r1*sin(theta); x_circ2 = r2*cos(theta); y_circ2 = r2*sin(theta);
end