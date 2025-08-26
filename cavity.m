% --- SETUP --- %

% Constants
c = 3*10^8; % speed of light, [m/s]
eps0 = (1/(36*pi))*10^(-9); % vacuum permittivity, [F/m]

% Adjustable parameters
L = 1; % length of cavity, [m]
D1 = 0.0254; % diameter of mirror 1, [m]
D2 = D1; % diameter of mirror 2, [m]
Rc1 = 3; % radius of curvature for mirror 1, [m]
Rc2 = Rc1; % radius of curvature for mirror 2, [m]
N = 1000; % number of mesh points along each dim of mesh grid
lambda = 780e-9; % laser wavelength, [m]
W = 5e-2; % domain half width, [m]
CFL = 0.0625;

% Grid
x = linspace(-W,W,N);
y = x;
dx = x(2) - x(1);
dy = dx;
%dz = CFL*4*k0*dx^2; % CFL-like condition, [m]
dz = L; % make each step a trip across the cavity
[X,Y] = meshgrid(x,y); % space domain

% Derived parameters
k0 = 2*pi/lambda; % freespace wavenumber, [m^-1]
D = 1i/(2*k0); % diffraction operator
mmask_R = (X.^2 + Y.^2 <= (D1/2)^2); % mirror 1 mask (RHS)
mmask_L = (X.^2 + Y.^2 <= (D2/2)^2); % mirror 2 mask (LHS)

% Set up mirror physical parameters for plotting
r1 = D1/2; % radius of mirror 1
r2 = D2/2; % radius of mirror 2
theta = linspace(0,2*pi,400);
x_circ1 = r1*cos(theta);
y_circ1 = r1*sin(theta);
x_circ2 = r2*cos(theta);
y_circ2 = r2*sin(theta);

% Input beam
w0 = 0.01; % input beam radius, [m]
E0 = exp(-(X.^2+Y.^2)/w0.^2); % input wave
E = E0;
I0 = 0.5*eps0*c*abs(E0).^2;

% Stability
g1 = 1 - L/Rc1; % stability parameter 1
g2 = 1 - L/Rc2; % stability paramter 2
g = g1*g2; % stability product, 0 < g < 1 for a stable cavity

% Set up frequency space
kx = (2*pi/(N*dx)) * (-N/2 : N/2-1);   % range from -pi/dx to +pi/dx
ky = kx;  % symmetric, since dx = dy
[KX, KY] = meshgrid(kx, ky);

% Free space transfer function of propagation 
H = exp(1i/(2*k0)*dz*(KX.^2+KY.^2));

% Mirror phase screen
M = exp(-1i*k0*(X.^2+Y.^2)/(Rc1));
%M = exp(1i*k0*(X.^2+Y.^2)/(2*Rc1));

%%
% --- SIMULATION --- %

fig = figure; 
subplot(1,2,1)
s=surf(X, Y, I0, 'LineStyle','none', 'DisplayName', '_nolegend_');
set(get(get(s,'Annotation'),'LegendInformation'), 'IconDisplayStyle','off');
hold on;
plot3(x_circ1, y_circ1, zeros(size(x_circ1)), 'r-', 'LineWidth', 2);
hold off;
title('Laser Mode at z = 0 m')
zlabel('Intensity')
xlabel('X [m]')
ylabel('Y [m]')

% Beam interacts with RHS mirror at L/2
E = E.*M;

save_interval = 1; % save frequency

% Propagation from L/2 (RHS) to -L/2 (LHS)
n1=0; % initialize step
l = L/2; % set current position
for z=0:dz:L % take steps of size dz, from 0 to L
    n1=n1+1; 
    Zp(n1) = z+dz; % total distance propagated so far
    l = l-dz; % current position within the cavity
    % Propagation in frequency domain, step 2
    % Gaussian Beam in frequency domain 
    FE = fft2(E); 
    % Propagated Gaussian beam in frequency domain 
    FE = FE.*fftshift(H); 
    % Propagated Gaussian beam in space domain 
    E = ifft2(FE); 
    % Step propagation through medium 
    if mod(n1, save_interval) == 0
        step_label = sprintf('step_%d', n1);
        Es.(step_label) = E; % save intermediate field
    end
end

I = 0.5*eps0*c*abs(E).^2;
subplot(1,2,2)
s=surf(X, Y, I, 'LineStyle','none');
set(get(get(s,'Annotation'),'LegendInformation'), 'IconDisplayStyle','off');
hold on;
plot3(x_circ2, y_circ2, zeros(size(x_circ2)), 'r-', 'LineWidth', 2, 'DisplayName', 'Mirror 2');
hold off;
title(sprintf('Laser Mode at z = %.1f m, After First Pass', L))
zlabel('Intensity')
xlabel('X [m]')
ylabel('Y [m]')

% Pass 2

% Beam clipped by LHS mirror at -L/2
E = E.*mmask_L;

% Beam interacts with LHS mirror at -L/2
E = E.*M;

fig = figure;
I = 0.5*eps0*c*abs(E).^2;
subplot(1,2,1)
surf(X, Y, I, 'LineStyle','none');
hold on;
plot3(x_circ2, y_circ2, zeros(size(x_circ2)), 'r-', 'LineWidth', 2, 'DisplayName', 'Mirror 2');
hold off;
title(sprintf('Laser Mode at z = %.1f m, Before Second Pass', L))
zlabel('Intensity')
xlabel('X [m]')
ylabel('Y [m]')

% Propagation from -L/2 back to L/2
n1=0; % initialize step
l = -L/2; % set current position
for z=0:dz:L % take steps of size dz, from 0 to L
    n1=n1+1; 
    Zp(n1) = z+dz; % total distance propagated so far
    l = l+dz; % current position within the cavity
    % Propagation in frequency domain, step 2
    % Gaussian Beam in frequency domain 
    FE = fft2(E); 
    % Propagated Gaussian beam in frequency domain 
    FE = FE.*fftshift(H); 
    % Propagated Gaussian beam in space domain 
    E = ifft2(FE); 
    % Step propagation through medium 
    if mod(n1, save_interval) == 0
        step_label = sprintf('step_%d', n1);
        Es.(step_label) = E; % save intermediate field
    end
end

% Beam clipped by RHS mirror at L/2
E = E.*mmask_R;

I = 0.5*eps0*c*abs(E).^2;
subplot(1,2,2)
surf(X, Y, I, 'LineStyle','none');
hold on;
plot3(x_circ2, y_circ2, zeros(size(x_circ2)), 'r-', 'LineWidth', 2, 'DisplayName', 'Mirror 2');
hold off;
title(sprintf('Laser Mode at z = %.1f m, After Second Pass', 2*L))
zlabel('Intensity')
xlabel('X [m]')
ylabel('Y [m]')


%% Post-processing

% Calculations

I = 0.5*eps0*c*abs(E0).^2; % intensity, [W/m^2]

% Graphing

figure(3);
surf(X, Y, I, 'LineStyle','none');
title('Laser Mode')
zlabel('Intensity')
xlabel('X [m]')
ylabel('Y [m]')

%% Helper functions

% Fourier transforms
function G = ft(g, delta)
    G = fftshift(fft(fftshift(g))) * delta;
end

function g = ift(G, delta_f)
    g = ifftshift(ifft(ifftshift(G))) * length(G) * delta_f;
end


%% Example, Poon and Kim

%Simulation of Gaussian Beam Focused by a Lens Using BPM 
%Paramters suggested for simulation : 
% Ld (light wavelength) =0.633, wo (waist)=10, 
% dz(sample distance along z)=800, Z(total final distance away from lens) =40000, 
% f(focal length)= 16000

%Gaussian Beam 
N=255; % sampling number 
L=50*10^(-3); %Display area 
%Ld=input('wavelength of light in [micrometers] = ?'); 
Ld = 0.633;
Ld=Ld*10^(-6); 
ko=(2*pi)/Ld; 
%wo=input('Waist of Gaussian Beam in [mm] = ?'); 
wo = 10;
wo=wo*10^(-3); 
%dz=input('step size of z (dz) in [mm] = ?'); 
dz =800;
dz=dz*10^(-3); 
%Z=input('destination of z in [mm] = ? '); 
Z=40000;
%Z=Z*10^(-3); %Focal length of Lens 
%f=input('Focal length of lens in [mm]= ?'); 
f= 16000;
f=f*10^(-3);

% dx : step size 
dx=L/N; 
for n= 1:256 
    for m=1:256 %Space axis 
        x(m)=(m-1)*dx-L/2; 
        y(n)=(n-1)*dx-L/2; 
        %Frequency axis 
        Kx(m)=(2*pi*(m-1))/(N*dx)-((2*pi*(N-1))/(N*dx))/2; 
        Ky(n)=(2*pi*(n-1))/(N*dx)-((2*pi*(N-1))/(N*dx))/2; 
    end 
end

[X,Y]=meshgrid(x,y); 
[KX, KY]=meshgrid(Kx,Ky); 

%Gaussian Beam in space domain 
Gau_ini=(1/(wo*pi^0.5))*exp(-(X.^2+Y.^2)./(wo^2));

%Energy of the initial Gaussian beam
Energy_ini=dx*dx*sum(sum(abs(Gau_ini).^2))

%Lens Equation 
L=exp(1i*ko/(2*f)*(X^2+Y^2)); 

%Gaussian Beam passed through the lens 
Gau_ini=Gau_ini.*L; % multiplicative phase operator in real space, step 1

%Free space transfer function of propagation 
H=exp(1i/(2*ko)*dz*(KX.^2+KY.^2));

%Iterative Loop 
Gau=Gau_ini; 
n1=0; 
for z=0:dz:Z 
    n1=n1+1; 
    Zp(n1)=z+dz; 
    % Propagation in frequency domain, step 2
    %Gaussian Beam in Frequency domain 
    FGau=fft2(Gau); 
    %Propagated Gaussian beam in Frequency domain 
    FGau=FGau.*fftshift(H); 
    %Propagated Gaussian beam in space domain 
    Gau=ifft2(FGau); 
    %Step propagation through medium 
    Gau_pro(:,n1)=Gau(:,127); 
end

%Energy of the final propagated Gaussian beam, 
%to check conservation of energy 
Energy_pro=dx*dx*sum(sum(abs(Gau).^2))

%axis in mm scale 
x=x*10^3; 
y=y*10^3; 
Zp=Zp*10^3; 
MAXl=max(max(abs(Gau_ini))); 
MAX2=max(max(abs(Gau))); 
MAX=max([MAXl MAX2]);

figure(1); 
mesh(x,y,abs(Gau_ini)) 
title('Initial Gaussian Beam') 
xlabel('x [mm]') 
ylabel('y [mm]') 
axis([min(x) max(x) min(y) max(y) 0 MAX]) 
axis square

figure(2); 
mesh(x,y,abs(Gau)) 
title('Propagated Gaussian Beam at Z') 
xlabel('x [mm]') 
ylabel('y [mm]') 
axis([min(x) max(x) min(y) max(y) 0 MAX]) 
axis square

figure(3); 
for l=1:n1
plot3(x',Zp(l)*ones(size(x')),abs(Gau_pro(:,l))) 
hold on 
end 
axis([min(x) max(x) min(Zp) max(Zp)]) 
grid on 
title('Beam profile along z') 
xlabel('x [mm]') 
ylabel('z [mm]') 
hold off

figure(4); 
A=max(abs(Gau_pro));  
N_Gau_pro=abs(Gau_pro) ./ A; 
contour(Zp,x,N_Gau_pro, [ exp(-1) exp(-1)], 'k') 
grid on 
title('Beam waist along z') 
xlabel('z [mm]') 
ylabel('x [mm]')