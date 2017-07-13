% Dino de Raad; June 29, 2017
% This script generates the fields for arbitary transducer arrangements
clear all;
close all;
clc;

% the rayleigh integral
% p(x,y,z,t) = ( j*k*rho_0*c_0*u_0*exp(j*omega*t) / 2*pi ) * sum(exp(-j*k*R)/R)*dS
%
% p is pressure
% x,y,z are space coordinates
% t is time
% j is sqrt(-1)
% k is the wavenumber, 2*pi/lambda = omega/c
% rho_0 is density
% c_0 is speed of sound
% u_0 is displacement amplitude
% omega is freqency (angular), or 2*pi*f
% R is the euclidean space between source and field points
% R = sqrt((x-x')^2 + (y-y')^2 + z^2), where x' and y' are source points

%% constants
f = 3e6; % in Hz

f = 1.057e6;

w = 2*pi*f; % in radians
c = 1500; % m/s in water @~25 deg C
k = w/c; 
rho = 1000; % kg/m^3 
u = 5.3e-5; % meters


%% source

R = 1e-2; % radius of curvature (m)
radius = 1e-3; % radius of aperture (m)

R = 44.5e-3; 
radius = 44.5e-3 / 2;

%% the Rayleigh length
% 
% RL = S/lambda
% S is the moving area of the source
% lambda is the wavelength

lambda = c/f; % in meters
S = 2 * pi*R^2 * (1 + sqrt(1-(radius./R).^2)); % in meters^2
RL = S/lambda; % in meters

%% the region S
% S is the surface of the vibrating piston
% the surface is a hemispherical shell of radius R with
% open face cut at the aperture radius

% number proportional to complexity of sphere (detail below)
% means higher numbers are required for smaller apertures
v_n = 200;
[xs, ys, zs] = sphere(v_n); % efficiency of this is ~O(n^2)
% can use function ellipsoid when geometry becomes more complicated

% Apply curvature radius
xs = xs(:) * R;
ys = ys(:) * R;
zs = zs(:) * R;

% Bowl
mask1 = (zs < 0);
xs = xs(mask1);
ys = ys(mask1);
zs = zs(mask1);

% Aperture
mask2 = (xs.^2 + ys.^2 <= radius^2);
xs = xs(mask2);
ys = ys(mask2);
zs = zs(mask2);

% Cuts
% Polynomials above the hemisphere
% mask3 = 

%% setup
theta = 0; % angle about x axis. 
phi = 0;   % angle about y axis
gamma = 0; % angle about z axis
Rx = @(t)[1 0 0; 0 cos(t) -sin(t); 0 sin(t) cos(t)]; % x rotation matrix for angle t
Ry = @(p)[cos(p) 0 sin(p); 0 1 0; -sin(p) 0 cos(p)]; % y ''
Rz = @(g)[cos(g) -sin(g) 0; sin(g) cos(g) 0; 0 0 1]; % z ''

num_z = 1;
maximum_angle = 2*pi;

XYZ = [xs.'; ys.'; zs.'];

if phi > 2 * pi || phi < 0
    error('Phi should be between 0 and 2pi')
end

if phi ~= 0
    XYZ = Ry(phi) * XYZ;
    
    if num_z > 1
        for l = 1:num_z
            XYZ = [XYZ Rz(l*maximum_angle/num_z)*XYZ];
        end
    end
end

%{
% Display S
scatter3(xs, ys, zs)

%axis([-radius radius -radius radius -2*R R])
figure
%}

%
%% field points

dzf = lambda / 10 * 2;
dxf = dzf * radius/R;
delta = .5;
zf = [-delta*R:dzf:-dzf 0:dzf:R*delta];
xf = [delta*-radius:dxf:-dxf 0:dxf:delta*radius];
yf = xf;
yf=0;

%xfyfzfR = Rx(theta)*[xf;yf;zf;];
%[xf,yf,zf] = squeeze(xfyfzfR(:,1,:));

% xy at z
% xz at y
% volume of points

%% Alternate Method
%
A = sqrt(-1) .* w .* rho .* c .* u/ (2 .* pi);
deltaS = dxf .^ 2;
R = @(x, y, z)sqrt((x-xs).^2 + (y-ys).^2 + (z-zs).^2);
pXYZ = @(x, y, z)sum(exp(-sqrt(-1).*k.*R(x, y, z))./R(x, y, z));

pf = zeros([length(xf), length(zf)]);

for n = 1:length(xf)
    for m = 1:length(zf)
        pf(m,n) = pXYZ(xf(n), 0, zf(m));
    end
    if mod(n, 10) == 0
        sprintf('%.0f %s', 100*n / length(xf), '%')
    end
end

pf = A .* deltaS .* pf;

%subplot(3,1,1)
waterfall(1000*xf,1000*zf,abs(pf)/1e3)
xlabel('x (mm)')
ylabel('z (mm)')
zlabel('P (Kp)')
hold on 
scatter3(1000*xs, 1000*zs, 1000*ys)
title('Pressure Field, Focus at 0')
%subplot(3,1,2)
%imagesc(1000*zf,1000*xf,real(pf))
%subplot(3,1,3)
%imagesc(1000*zf,1000*xf,imag(pf))
figure
subplot(2,1,1)
pCAf = pf(round(length(pf)/2), :);
pAf = pf(:, round(length(pf)/2));
plot(1000*xf, real(pAf))
hold on
plot(1000*xf, imag(pAf))
hold on
plot(1000*xf, abs(pAf), 'k')
[~, I] = max(abs(pCAf));
title(sprintf('Axial Pressure | Max: %f mm', 1000*zf(I)))
subplot(2,1,2)
plot(1000*zf, real(pCAf))
hold on
plot(1000*zf, imag(pCAf))
hold on
plot(1000*zf, abs(pCAf), 'k')
[~, I] = max(abs(pCAf));
title(sprintf('Cross-Axial Pressure | Max: %f mm', 1000*xf(I)))
