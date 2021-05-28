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
   
    postParticles.x_r = ... % 1xN_particles matrix
    postParticles.y_r = ... % 1xN_particles matrix
    postParticles.phi = ... % 1xN_particles matrix
    postParticles.kappa = ... % 1xN_particles matrix
    
    % and leave the function
    return;
end % end init

%% Mode 2: Estimator iteration.
% If km > 0, we perform a regular update of the estimator.

% Implement your estimator here!


% Prior Update:


% Posterior Update:

postParticles.x_r = ...
postParticles.y_r = ...
postParticles.phi = ...
postParticles.kappa = ...

end % end estimator

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