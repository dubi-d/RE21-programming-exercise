simConst = EstimatorConst();
[basePts, vecs] = initContourVectors(simConst.contour);
% n_parts = 20000;
% tic
% [tmin] = calcMinDistances(basePts, vecs, ones(n_parts,1)*[2], ones(n_parts,1)*[0.5], ones(n_parts,1)*[5*pi/4], -0.2:0.4/(n_parts):0.2);
% toc

n_particles = 20;
px = [rand(1, n_particles)*3, 2, 1.5];
py = [rand(1, n_particles)*3, 2, 1.5];
phi = [rand(1, n_particles)*2*pi, pi, pi];
kappa = [zeros(1,n_particles), 0.1, -0.1];
tmin = calcMinDistances(basePts, vecs, px, py, phi, kappa);

% Plot
figure(1)
clf
% Plot contour
spacing = 0.6;
contour = simConst.contour;
contour(8, 1) = contour(8, 1);
contour(9, 1) = contour(9, 1);
axis([min(contour(:,1))-spacing,max(contour(:,1))+spacing,...
      min(contour(:,2))-spacing,max(contour(:,2))+spacing])
hold on
for i = 1:size(contour,1)-1
    plot([contour(i,1),contour(i+1,1)],[contour(i,2),contour(i+1,2)],'k--')
end
plot([contour(end,1),contour(1,1)],[contour(end,2),contour(1,2)],'k--')

% Plot particles
scatter(px, py, 'bo')
for i = 1:length(px)
    line([px(i); px(i)+10*cos(phi(i))], [py(i); py(i)+10*sin(phi(i))])
    scatter(px(i)+tmin(i).*cos(phi(i)), py(i)+tmin(i).*sin(phi(i)), 'rx')
end

hold off
axis([min(contour(:,1))-spacing,max(contour(:,1))+spacing,...
      min(contour(:,2))-spacing,max(contour(:,2))+spacing])

%%

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
