%%
%%
%%
%% Code to test tilt and shift operations observed in rotating
%% frame of reference 
%% CL 10/18/2025

clear all;
clc;
figure;

v = VideoWriter('cavitymode_Om0.005.avi','Motion JPEG AVI');
v.Quality = 95;
open(v);
 
 Lx=2; %Length of square transverse domain (one side) [m]
 Ld=1064*1e-9; %Laser wavelength [m]
 wo=1e-1; %initial beam width [m]
 ko=(2*pi)/Ld; %wavenumber 
 Zmax=20e3; %destination of z in [m] 
 Z0 = -Zmax; %Starting location
 L = Zmax - Z0;
 Nz = 1e2;
 dz = L/Nz;
 Omega = -0.01; %rad/sec
 c = 3e8; %m/s

%N : sampling number
N=511;  
dx=Lx/N; % dx : step size 

%Create grid
for n=1:N+1 
    for m=1:N+1
%Space axis
x(m)=(m-1)*dx-Lx/2;
y(n)=(n-1)*dx-Lx/2;
    end
end
[X,Y]=meshgrid(x,y);


%Gaussian Beam in space domain
Gau_ini=(1/(wo*pi*0.5))*exp(-(X.^2+Y.^2)./(wo^2));

%Beam Propagation in Inhomogeneous
%Iterative Loop 
Gau=Gau_ini; 
%Centroid locations to begin
centerx(1) = trapz(trapz(X.*abs(Gau).^2))/trapz(trapz(abs(Gau).^2));
centery(1) = trapz(trapz(Y.*abs(Gau).^2))/trapz(trapz(abs(Gau).^2));
%z locations to evaluate field. 
z = linspace(Z0,Zmax,Nz);

for i = 1:Nz-1

%Propagation using angular spectrum method (Schmidt et al.)
[x2dum, Gau] = ang_spec_prop(Gau,Ld,dx,dx,dz);

%new X-position after trajectory along characteristic curve
Xnew = X - 0.5*Omega/c * (z(i+1).^2-z(i).^2);

%Analytic integration of phase term along characteristic 
Xint = X*(z(i+1)-z(i)) - 0.5*Omega/c * ...
    (1/3 * (z(i+1).^3-z(i)^3) + 0.5*z(i)*(z(i+1)-z(i)));

%Phase screen implementation for tilt term
A = exp(-1i * ko * Omega / c * Xint);
Gau = Gau .* A; %apply tilt

%Implementation of shift interpolation
Gau = interp2(Xnew,Y,Gau,X,Y,'spline');

centerx(i+1) = trapz(trapz(X.*abs(Gau).^2))/trapz(trapz(abs(Gau).^2)); 
centery(i+1) = trapz(trapz(Y.*abs(Gau).^2))/trapz(trapz(abs(Gau).^2)); 

imagesc(x,y,abs(Gau)); axis([-0.5 0.5 -0.5 0.5]); axis square; xlabel('x [m]'); ylabel('y [m]'); 
hold on; plot(centerx(i+1),centery(i+1),'ro'); hold off; frame = getframe(gcf); display(z(i));
writeVideo(v,frame);

end

x_ana = -Omega/c * (z.^2 + 0.5*L*z); %analytic solution from ray optics

%imagesc(abs(Gau)); title('Beam Profile'); getframe(); display(z(i));
%Everything below: Plotting


figure;
plot(centerx,z,'b',x_ana,z,'r'); xlabel('Distance [m]'); ylabel('Displacement [m]');


close(v);


