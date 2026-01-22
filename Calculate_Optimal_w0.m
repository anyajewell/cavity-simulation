% What mirror curvatures result in the smallest mode size on the mirror, as
% a function of the separation L?

% Settings
%w0 = 0.1; % initial beam waist, [m]
lambda = 1064e-9; % wavelength, [nm]
%L = 100; % cavity length, [km]
%L = L*1e3;
L = 10000:10:300000

% Calculations
w0 = sqrt(lambda*L./(2*pi)); % input waist that will give the smallest spot size on mirror
z = L./2; % looking at the spot size at the mirror (either one)
zR = (pi*w0.^2)./lambda; % Rayleigh range
wz = w0.*sqrt(1 + (z./zR).^2); % spot size at L/2
Rc = L./2 + zR.^2./(L./2) % radius of curvature required to minimize spot size

% Plotting
figure;
plot(L/1000, Rc./L, 'LineWidth', 1)
title({'Radius of Curvature Required to Minimize'}, {'Beam Spot Size on Mirrors'}, 'Interpreter','latex')
ylabel('$R_c / L$', 'Interpreter', 'latex');
xlabel('$L$ [km]', 'Interpreter', 'latex')

