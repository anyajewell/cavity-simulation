%------------Ray optics bouncing simulation----------%
%------------In a rotating cavity--------------------%
%------------10/13/25 C. Limbach---------------------%

clear;
hold on;

%% -------------Fundamental and Cavity Parameters -------------%% 
c = 3e8; %m/s
Omega = 0.0005; %0.001; %0.001;  %rad/sec rotation rate
L = 100e3; %m %cavity length

%Mirror radii of curvature
Rc1 = 1.5*L; %L; %Rc = 2L gives a stable config
Rc2 = 1.5*L;
tilt1 = 0; %1e-6; %rad/m  %-2 and 0.667 works
tilt2 = 0; %1e-6; %rad/m

%Derived parameters
g1 = 1-L/Rc1;
g2 = 1-L/Rc2;
display(g1*g2);

%% -----------------Initial conditions--------------------%%
th0 = 0; %-L*Omega/(2*c); %% 
% look-ahead angle is equal to -Omega*L/(2*c);
x0 = 0.01; %0.001
z0 = -L/2;
z = linspace(-L/2,L/2,100);
Nroundtrips = 100;

%% --------------------Main Loop----------------------------%%

for k = 1:Nroundtrips

%-------------Foward propagation----------------%

%Equation to implement
%𝑥−𝑥_0=Ω/𝑐 (𝑧^2−𝑧_0^2 )+(𝜃_0−(Ω𝑧_0)/𝑐)(𝑧−𝑧_0)
%th = th0 + Omega/c * (z(end)-z0)

x = x0 + Omega/c * (z.^2 - z0^2) + (th0 - Omega*z0/c)*(z-z0);
th0 = th0 + Omega/c * (z(end)-z0);

plot(x,z); getframe(); hold on; %pause(1); 

%----------------Reverse propagation------------------%
%Reset initial conditions
x0 = x(end);
th0 = th0 - 2*x0/Rc2; %mirror term
th0 = th0 - 2*tilt2;
z0 = L/2;

%Equation to implement:
%−(𝑥−𝑥_0)=Ω/𝑐 (𝑧^2−𝑧_0^2 )+(𝜃_0−(Ω𝑧_0)/𝑐)(𝑧−𝑧_0)
%th = th0 + Omega/c * (z(1)-z0)

x = x0 - Omega/c * (z.^2 - z0^2) - (th0-(Omega*z0/c))*(z-z0);
th0 = th0 + Omega/c * (z(1)-z0);

plot(x,z); getframe(); %pause(1);

%%     ------------------------------------------- %%
%Reset initial conditions
x0 = x(1);
th0 = th0 - 2*x0/Rc1; %Angle change from Mirror radius of curvature of mirror
th0 = th0 - 2*tilt1;
z0 = -L/2;

%pause(1); %slow down simulation

end

