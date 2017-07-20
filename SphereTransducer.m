% Generates a simple spherical trasducer.
% rK is the radius of curvature (m), rA is the radius of aperture (m) 
% opional f (Hz) to increase resolution ([] defaults), optional 
% rH (m) to put a hole in the center.
function S = SphereTransducer(rK, rA, varargin)
    % the surface is a hemispherical shell of radius R with
    % open face cut at the aperture radius
    
    % number proportional to complexity of sphere (detail below)
    % means higher numbers are required for smaller apertures
    v_n = 300;
    rH = 0;
    
    switch length(varargin)
        case 0
        case 1
            if ~isempty(varargin{1})
                v_n = ceil(300/1e6 * varargin{1});
            end
        case 2
            if ~isempty(varargin{1})
                v_n = ceil(300/1e6 * varargin{1});
            end
            rH = varargin{2};
        otherwise
            error('Too many inputs')
    end

    
    [xs, ys, zs] = sphere(v_n); % efficiency of this is ~O(n^2)
    % can use function ellipsoid when geometry becomes more complicated

    % Apply curvature radius
    xs = xs(:) * rK;
    ys = ys(:) * rK;
    zs = zs(:) * rK;

    % Bowl
    mask1 = (zs < 0);
    xs = xs(mask1);
    ys = ys(mask1);
    zs = zs(mask1);

    % Aperture
    mask2 = (xs.^2 + ys.^2 <= rA^2);
    xs = xs(mask2);
    ys = ys(mask2);
    zs = zs(mask2);
    
    % Hole
    mask3 = (xs.^2 + ys.^2 >= rH^2);
    xs = xs(mask3);
    ys = ys(mask3);
    zs = zs(mask3);
    
    S = [xs ys zs];
end