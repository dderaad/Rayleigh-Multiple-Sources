%% Demo 2

close all; clc; clear all;

r = 44.5e-3;
S = SphereTransducer(r, r/2, [], r/4);
x_lim = [-r/2 r/2];
z_lim = [-r/2 r/2];
f = 1.057e6;
res = log(10)/log(10);

S_list = cell([4 1]);
S_list{1} = S;


i = 2;
for g = [-pi/4 pi/4]
    S_list{i} = rotate(S,-pi/4,g,0, [1 2 3]);
    i = i + 1;
end

for g = [0]
    S_list{i} = rotate(S,-7*pi/12,g,0, [2 1 3]);
    i = i + 1;
end

xS = [];
yS = [];
zS = [];

for n = 1:length(S_list)
    S_temp = S_list{n};
    xS = [xS; S_temp(:,1)];
    yS = [yS; S_temp(:,2)];
    zS = [zS; S_temp(:,3)];
end


C = 127*(zS - min(zS))/(max(zS) - min(zS));
scatter3(1000*xS, 1000*zS, 1000*yS, 25, C);
colormap(hot)
set(gca,'Color',[0.925 0.925 0.925]);
%{
hold on
scatter3(1000*S(:,1), 1000*S(:,3), 1000*S(:,2))
%}

axis image
drawnow;
hold on

P = 0;

for n = 1:length(S_list)
    [p, x, z] = RI(S_list{n}, x_lim, z_lim, res, f);
    P = P + abs(p);
end

sf = 1000 * sqrt(sum(range([xS yS zS]).^2))/2;
surf(1000*x, 1000*z, sf*P/max(max(P)))
title(sprintf('Gain: %.2g Pa', max(max(P))))
shading interp;