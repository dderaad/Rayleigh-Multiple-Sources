%% Demo 1

close all; clc;
r = 44.5e-3;
S = SphereTransducer(r, r/2);
x_lim = [-r*1.125*.5 r*1.125*.5];
z_lim = [-r/2 r/2];
f = 1.057e6;

%f = 1.057e6
%x_lim = [-r r]; 
%z_lim = [-r*.8 r*1.2];

[p, x, z] = RI(S, x_lim, z_lim, .7, f);
P = abs(p);

%% Real dimensions of transducer to focus.
figure('Color',[0.9 0.9 0.9])

xS = S(:,1);
yS = S(:,2);
zS = S(:,3);
C = 10*(zS - min(zS))/(max(zS) - min(zS));
scatter3(1000*xS, 1000*zS, 1000*yS, 25, C, 'filled')
colormap(hot(10))
hold on 
sf = 1000 * sqrt(sum(range(S).^2))/2;
surf(1000*x,1000*z,sf*P/max(max(P)))
xlabel('x (mm)')
ylabel('z (mm)')
zlabel('y (mm), P (arb.)')
axis equal
set(gca,'ztick',[])
set(gca,'zticklabel',[])
title('Pressure Field and Transducer, Focus at 0')
shading interp;
set(gca,'Color',[0.925 0.925 0.925]);

%% Axial, Cross-Axial Response
figure('Color',[0.9 0.9 0.9])

subplot(2,1,1)
[I, J] = find(P==max(max(P)));
pCAf = P(I, :);
pAf  = P(:, J);
plot(1000*z, abs(pAf), 'k')
xlabel('axis (mm)')
ylabel('P (Pa)')
[M_A, I_A] = max(abs(pAf));
title(sprintf('Axial Pressure | Max: %.2g (Pa) @ %.1f mm', M_A, 1000*z(I_A)))
xlim([1000*z(1) 1000*z(end)])
ylim([0, M_A * 1.125])

subplot(2,1,2)
plot(1000*x, abs(pCAf), 'k')
xlabel('cross-axis (mm)')
ylabel('P (Pa)')
[M_CA, I_CA] = max(abs(pCAf));
title(sprintf('Cross-Axial Pressure | Max: %.2g (Pa) @ %.1f mm',...
	M_CA, 1000*x(I_CA)))
ylim([0, M_CA * 1.125])
xlim([1000*x(1) 1000*x(end)])
suptitle(sprintf('Overall Gain: %.2g Pa', max([M_A M_CA])))

%% Contour of focal region.
figure('Color',[0.9 0.9 0.9])

[C, h] = contour(1000*x,1000*z,P,[M_A/2 M_CA/2]);
C = C(:,2:end);
xC = C(1,:);
zC = C(2,:);
plot(xC, zC)
hold on
xh = [-max(xC)*1.125 max(xC)*1.125]; yh = [1000*z(I_A) 1000*z(I_A)];
xv = [1000*x(I_CA) 1000*x(I_CA)]; yv = [-max(zC)*1.125 max(zC)*1.125];
plot(xh,yh,xv,yv,'k')

xlim([-max(xC)*1.125 max(xC)*1.125])
ylim([-max(zC)*1.125 max(zC)*1.125])
xlabel('cross-axis (mm)')
ylabel('axis (mm)')

% determine the axial and cross-axial distances
% A
yvC = interp1(yv, xv, xC, 'spline');
idx = find(abs(yvC - xC) < .1);
%hold on
%scatter(xInter, yInter)

dA = sqrt((min(xInter) - max(xInter))^2+(min(yInter) - max(yInter))^2);
% CA
xvC = interp1(xh, yh, zC, 'spline');
idx = find(abs(xvC - zC) < .1);
xInter = xC(idx);
yInter = zC(idx);

dCA = sqrt((min(xInter) - max(xInter))^2+(min(yInter) - max(yInter))^2);
fA = polyarea(xC, zC); % Area of focal region.

axis image
title(sprintf(...
	'Focal Area: %.3g mm^2 | Axial Distance: %.3g mm | Cross-Axial Distance: %.3g mm'...
	, fA, dA, dCA))

fV = pi*trapz(zC(xC>-eps), xC(xC>-eps).^2); % Area of revolution.

suptitle(sprintf('Focal Volume (Axisymmetric Source): %.3g mm^3', fV))
