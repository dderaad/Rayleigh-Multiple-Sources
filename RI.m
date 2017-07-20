% Computes the Rayleigh integral for the given transducer (pre-focused at 0)
% in a square specified by x_lim and z_lim. Res of 1 computes the field
% at a resolution of ten times per wavelength. Evey integer increment of Res
% will increase the field resolution tenfold.
% Assumes a water medium @ 25 deg C and 1MHz frq (if not specified).
% Optional arguments: f (frequency, Hz), rho (density, kg/m^3),
% c (sound speed, m/s)
function [field, xf, zf] = RI(S, x_lim, z_lim, Res, varargin)
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

    f = 1e6; % in Hz
    c = 1500; % m/s in water @~25 deg C
    rho = 1000; % kg/m^3

    switch length(varargin)
        case 0
        case 1
            f = varargin{1};
        case 2
            f = varargin{1};
            c = varargin{2};
        case 3
            f = varargin{1};
            c = varargin{2};
            rho = varargin{3};
        otherwise
            error('Too many inputs.')
    end
    
    
    w = 2*pi*f; % in radians / s
    k = w/c;  
    lambda = c/f;
    
    sizeS = size(S);
    if sizeS(1) == 3
        S = S.';
    elseif sizeS(2) == 3
    else
        error('Source points not in valid format: (X,Y,Z), (X,Y,Z).''')
    end
    
    xs = S(:,1);
    ys = S(:,2);
    zs = S(:,3);
    
    mins = min(S);
    maxs = max(S);
    radii = (maxs - mins)/2;
    r = min(radii);
    R = mean(radii(radii~=r));
    
    if R ~= 0
       [~, E] = ellipke(1-r/R);
       rA = R*E; % This is the half arclength of the sphere
    else
       rA = 1;
    end
    
    u = 2./(k*rho*c*rA); % m (initial particle velocity amplitude)
    
    dzf = lambda / 10^Res;
    dxf = dzf * (x_lim(2)-x_lim(1))/(z_lim(2)-z_lim(1));
    zf = [z_lim(1):dzf:-dzf 0:dzf:z_lim(end)];
    xf = [x_lim(1):dxf:-dxf 0:dxf:x_lim(end)];
    
    A = sqrt(-1) .* w  .* rho .* c .* u ./ (2 .* pi);
    deltaS = approximate_dS(S).^2;
    R = @(x, y, z)sqrt((x-xs).^2 + (y-ys).^2 + (z-zs).^2);
    
    field = zeros([length(xf), length(zf)]);
    
    for n = 1:length(xf) 
        Rn = R(xf(n), 0, zf);
        field(n, :) = A .* deltaS .* sum(exp(-sqrt(-1).*k.*Rn)./Rn);
    end
    field = field.';
end

% Works for any arbitrary arrangement
function ds = approximate_dS(S)
    sample = 1:50;
    ds = zeros(length(sample), 1);
    for n = sample
        I = datasample(1:length(S), ceil(.1*length(S)));
        points = S(I, :);
        potnn = S;
        potnn(I, :) = [];

        [~, D] = knnsearch(potnn, points);
        ds(n) = mean(D);
    end
    ds = median(ds);
end