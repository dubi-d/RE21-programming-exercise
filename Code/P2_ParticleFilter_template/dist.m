const = EstimatorConst();
[basePts, vecs] = initContourVectors(const.contour);
n_parts = 20000;
tic
% [tmin] = calcMinDistances(basePts, vecs, ones(n_parts,1)*[2], ones(n_parts,1)*[0.5], ones(n_parts,1)*[5*pi/4], -0.2:0.4/(n_parts):0.2);
tmins = calcMinDistances(basePts, vecs, [1], [1], [0], [0])
toc

function [basePts, vecs] = initContourVectors(contour)
    contour = [contour; contour(1,:)];
    vecs = diff(contour);
    basePts = contour(1:end-1, :);
end

function [basePts, vecs] = updateContour(basePts, vecs, kappa)
    basePts(8:9, 1) = kappa;
    vecs(7, 1) = vecs(7,1) + kappa;
    vecs(9, 1) = vecs(9,1) - kappa;
end

function [tmin] = calcMinDistParticle(basePts, vecs, px, py, phi, kappa)
    [basePts, vecs] = updateContour(basePts, vecs, kappa);
    
    p = [px, py];
    r = [cos(phi), sin(phi)];
    
    t = diff((p - basePts) .* fliplr(vecs), 1, 2) ./ diff(fliplr(r) .* vecs, 1, 2);
    s = diff((basePts - p) .* fliplr(r), 1, 2) ./ diff(fliplr(vecs) .* r, 1, 2);
    
    tmin = min(t((s >= 0) & (s <= 1) & (t >= 0)));
    if isempty(tmin) 
        tmin = inf;
    end
end

function t_mins = calcMinDistances(basePts, vecs, px, py, phi, kappa)
    n_pts = length(phi);
    t_mins = zeros(n_pts, 1);
    for i = 1:n_pts
        t_mins(i) = calcMinDistParticle(basePts, vecs, px(i), py(i), phi(i), kappa(i));
    end
end
