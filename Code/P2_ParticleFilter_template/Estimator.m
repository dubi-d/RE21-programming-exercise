function [postParticles] = Estimator(prevPostParticles, sens, act, estConst, km)
% The estimator function. The function will be called in two different
% modes: If km==0, the estimator is initialized. If km > 0, the
% estimator does an iteration for a single sample time interval using the 
% previous posterior particles passed to the estimator in
% prevPostParticles and the sensor measurement and control inputs.
%
% Inputs:
%   prevPostParticles   previous posterior particles at time step k-1
%                       The fields of the struct are [1xN_particles]-vector
%                       (N_particles is the number of particles) 
%                       corresponding to: 
%                       .x_r: x-locations of the robot [m]
%                       .y_r: y-locations of the robot [m]
%                       .phi: headings of the robot [rad]
%                       .kappa: wall offset [m]
%                           
%   sens                Sensor measurement z(k), scalar
%
%   act                 Control inputs u(k-1), [1x2]-vector
%                       act(1): u_f, forward control input
%                       act(2): u_phi, angular control input
%
%   estConst            estimator constants (as in EstimatorConst.m)
%
%   km                  time index k, scalar
%                       corresponds to continous time t = k*Ts
%                       If km==0 initialization, otherwise estimator
%                       iteration step.
%
% Outputs:
%   postParticles       Posterior particles at time step k
%                       The fields of the struct are [1xN_particles]-vector
%                       (N_particles is the number of particles) 
%                       corresponding to: 
%                       .x_r: x-locations of the robot [m]
%                       .y_r: y-locations of the robot [m]
%                       .phi: headings of the robot [rad]
%                       .kappa: wall offset [m]
%
%
% Class:
% Recursive Estimation
% Spring 2021
% Programming Exercise 2
%
% --
% ETH Zurich
% Institute for Dynamic Systems and Control
% Raffaello D'Andrea, Matthias Hofer, Carlo Sferrazza
% hofermat@ethz.ch
% csferrazza@ethz.ch

% Set number of particles:
N_particles = 10; % obviously, you will need more particles than 10.

%% Mode 1: Initialization
if (km == 0)
    % Do the initialization of your estimator here!
    
    [x_samples, y_samples] = sample_init_positions(N_particles, estConst);
    postParticles.x_r = x_samples; % 1xN_particles matrix
    postParticles.y_r = y_samples; % 1xN_particles matrix
    postParticles.phi = get_uni_vec(estConst.phi_0, N_particles); % 1xN_particles matrix
    postParticles.kappa = get_uni_vec(estConst.l, N_particles); % 1xN_particles matrix
    
    % and leave the function
    return;
end % end init

%% Mode 2: Estimator iteration.
% If km > 0, we perform a regular update of the estimator.

% Implement your estimator here!


% Prior Update:
[x_r, y_r, phi, kappa] = move_particles(prevPostParticles, act, estConst, N_particles);

% Posterior Update:
[basePts, vecs] = initContourVectors(estConst.contour);
t_min = calcMinDistances(basePts, vecs, x_r, y_r, phi, kappa);
p_meas = measurementProbabilities(t_min, sens, estConst.epsilon);

weights = p_meas/sum(p_meas);


% % postParticles.x_r = ...
% % postParticles.y_r = ...
% % postParticles.phi = ...
% % postParticles.kappa = ...

end % end estimator

function [uni] = get_uni_vec(bound, N)
% get a vector (1xN) of samples from uniform(-bound, bound)
% (shift and scale standard uniform distribution)
uni = (rand(1, N) - 0.5) * 2 * bound;
end

function [x_samples, y_samples] = sample_init_positions(N, estConst)
    pA = estConst.pA; % center of circle A
    pB = estConst.pB; % center of circle B
    r = estConst.d; % radius of circles

    x_samples = zeros(1, N);
    y_samples = zeros(1, N);

    for i = 1:N
        isA = (rand < 0.5); % choose circle
        xi = (rand - 0.5) * 2 * r; % sample x in generic circle
        yBound = sqrt(r^2 - xi^2); % compute range of posible y given x
        yi = (rand - 0.5) * 2 * yBound; % sample y
        x_samples(i) = xi + pA(1) * isA + pB(1) * (~isA); % shift to given circle
        y_samples(i) = yi + pA(2) * isA + pB(2) * (~isA);
    end
end

function [x_r, y_r, phi, kappa] = move_particles(part, act, estConst, N)
    % process noise is ~Uniform(-bound, bound)
    f_bound = estConst.sigma_f;
    phi_bound = estConst.sigma_phi;

    % vectorize inputs
    u_f = repmat(act(1), 1, N);
    u_phi = repmat(act(2), 1, N);

    % apply process model including noise
    x_r = part.x_r + (u_f + get_uni_vec(f_bound, N)) .* cos(part.phi);
    y_r = part.y_r + (u_f + get_uni_vec(f_bound, N)) .* sin(part.phi);
    phi = part.phi + u_phi + get_uni_vec(phi_bound, N);
    kappa = part.kappa;
end

function [particlesResampled, weightsResampled] = resample(particles, weights)

    % systematic resampling (low variance)

    % weights are already normalized beforehand (normalization of weights is required or this won't work)

    numSpokes = size(particles, 1);  % number of spokes of the resampling wheel (= number of particles to sample)
    u = rand() / numSpokes;  % first spoke's position along the arc of the wheel
    sumWeights = weights(1);  % initialize accumulating sum of weights (= arc length covered so far)
    j = 1;
    particlesResampled = zeros(size(particles));
    weightsResampled = zeros(size(weights));

    % going through all the spokes
    for i = 1:numSpokes

       % check which particle (= arc segment with length of its weight) gets hit by the current spoke
       while sumWeights < u
           j = j + 1;
           sumWeights = sumWeights + weights(j);

       end

       % select the sample that gets hit by the current spoke
       particlesResampled(i, :) = particles(j, :);
       weightsResampled(i) = weights(j);

       u = u + 1/numSpokes;  % move on to the next spoke

    end

    % re-normalize the new weights
    weightsResampled = weightsResampled / sum(weightsResampled);

end

%% Distance calculation
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
end

function t_mins = calcMinDistances(basePts, vecs, px, py, phi, kappa)
    n_pts = length(phi);
    t_mins = zeros(n_pts, 1);
    for i = 1:n_pts
        t_mins(i) = calcMinDistParticle(basePts, vecs, px(i), py(i), phi(i), kappa(i));
    end
end

%% Measurement Probability

function measProb = measurementProbabilities(t, z, epsilon)
    d = abs(t-z);
    
    measProb = zeros(size(d));
    measProb(d < 2*epsilon) = 1/(5*epsilon)*(2 - d(d < 2*epsilon)/epsilon);
    measProb(d > 2*epsilon & d < 3*epsilon) = 1/(5*epsilon)*(1- 2*abs(d(d > 2*epsilon & d < 3*epsilon)-2.5*epsilon)/epsilon);
end