function [posEst,linVelEst,oriEst,windEst,driftEst,...
          posVar,linVelVar,oriVar,windVar,driftVar,estState] = ...
    Estimator(estState,actuate,sense,tm,estConst)
% [posEst,linVelEst,oriEst,windEst,driftEst,...
%    posVar,linVelVar,oriVar,windVar,driftVar,estState] = 
% Estimator(estState,actuate,sense,tm,estestConst)
%
% The estimator.
%
% The function is initialized for tm == 0, otherwise the estimator does an
% iteration step (compute estimates for the time step k).
%
% Inputs:
%   estState        previous estimator state (time step k-1)
%                   May be defined by the user (for example as a struct).
%   actuate         control input u(k-1), [1x2]-vector
%                   actuate(1): u_t, thrust command
%                   actuate(2): u_r, rudder command
%   sense           sensor measurements z(k), [1x5]-vector, INF entry if no
%                   measurement
%                   sense(1): z_a, distance measurement a
%                   sense(2): z_b, distance measurement b
%                   sense(3): z_c, distance measurement c
%                   sense(4): z_g, gyro measurement
%                   sense(5): z_n, compass measurement
%   tm              time t_k, scalar
%                   If tm==0 initialization, otherwise estimator
%                   iteration step.
%   estestConst        estimator estConstants (as in EstimatorestConst.m)
%
% Outputs:
%   posEst          position estimate (time step k), [1x2]-vector
%                   posEst(1): p_x position estimate
%                   posEst(2): p_y position estimate
%   linVelEst       velocity estimate (time step k), [1x2]-vector
%                   linVelEst(1): s_x velocity estimate
%                   linVelEst(2): s_y velocity estimate
%   oriEst          orientation estimate (time step k), scalar
%   windEst         wind direction estimate (time step k), scalar
%   driftEst        estimate of the gyro drift b (time step k), scalar
%   posVar          variance of position estimate (time step k), [1x2]-vector
%                   posVar(1): x position variance
%                   posVar(2): y position variance
%   linVelVar       variance of velocity estimate (time step k), [1x2]-vector
%                   linVelVar(1): x velocity variance
%                   linVelVar(2): y velocity variance
%   oriVar          variance of orientation estimate (time step k), scalar
%   windVar         variance of wind direction estimate(time step k), scalar
%   driftVar        variance of gyro drift estimate (time step k), scalar
%   estState        current estimator state (time step k)
%                   Will be input to this function at the next call.
%
%
% Class:
% Recursive Estimation
% Spring 2021
% Programming Exercise 1
%
% --
% ETH Zurich
% Institute for Dynamic Systems and Control
% Raffaello D'Andrea, Matthias Hofer, Carlo Sferrazza
% hofermat@ethz.ch
% csferrazza@ethz.ch

%% Initialization
if (tm == 0)
    % Do the initialization of your estimator here!
    
    % initial state mean
    posEst = zeros(1, 2); % 1x2 matrix
    linVelEst = zeros(1, 2); % 1x2 matrix
    oriEst = 0; % 1x1 matrix
    windEst = 0; % 1x1 matrix
    driftEst = 0; % 1x1 matrix
    
    % initial state variance
    posVar = [1/4*estConst.StartRadiusBound^2, 1/4*estConst.StartRadiusBound^2]; % 1x2 matrix
    linVelVar = zeros(1, 2); % 1x2 matrix
    oriVar = (2 * estConst.RotationStartBound)^2 / 12; % 1x1 matrix
    windVar = (2 * estConst.WindAngleStartBound)^2 / 12; % 1x1 matrix
    driftVar = (2 * estConst.GyroDriftStartBound)^2 / 12; % 1x1 matrix
    
    % estimator variance init (initial posterior variance)
    estState.Pm = diag([posVar, linVelVar, oriVar, windVar, driftVar]);
    % estimator state
    estState.xm = [posEst, linVelEst, oriEst, windEst, driftEst];
    % time of last update
    estState.tm = tm;
    return
end

%% Estimator iteration.
% get time since last estimator update
n_states = 7;
% dt = tm - estState.tm;

xm = estState.xm';
Pm = estState.Pm;

% prior update
% Intregrate states and variance simultaneously
xAndP_0 = [xm; reshape(Pm, [], 1)];
[~, xAndP] = ode45(@(t, xAndP) evalqAndPdot(xAndP, actuate, estConst, n_states), [estState.tm, tm], xAndP_0);
P_prior = reshape(xAndP(end, n_states+1:end), n_states, n_states);
x_prior = xAndP(end, 1:n_states)';


% measurement update
H = evalH(xm, estConst);
M = evalM();
R = diag([estConst.DistNoiseA, estConst.DistNoiseB, estConst.DistNoiseC,...
          estConst.GyroNoise, estConst.CompassNoise]);
K = (P_prior * H') / (H * P_prior * H' + M * R * M');

% No measurement when sense == inf -> remove this column from K
K(:,sense'==inf) = 0;
% Remove inf in sense as otherwise multiplication with 0 yields NaN
sense(sense==inf) = 0;

h = evalh(x_prior, estConst);
x_posterior = x_prior + K * (sense' - h);
P_posterior = (eye(n_states) - K*H) * P_prior * (eye(n_states) - K*H)' + K*M*R*M'*K';

% Update estState
estState.tm = tm; % update measurement update time
estState.xm = x_posterior';
estState.Pm = P_posterior;

% Get resulting estimates and variances
% Output quantities
posEst = x_posterior(1:2);
linVelEst = x_posterior(3:4);
oriEst = x_posterior(5);
windEst = x_posterior(6);
driftEst = x_posterior(7);
% 
posVar = diag(P_posterior(1:2, 1:2));
linVelVar = diag(P_posterior(3:4, 3:4));
oriVar = P_posterior(5, 5);
windVar = P_posterior(6, 6);
driftVar = P_posterior(7, 7);

end


function q = evalq(x, u, estConst)
    sx = x(3);
    sy = x(4); 
    phi = x(5);
    rho = x(6);
    Cdh = estConst.dragCoefficientHydr;
    Cda = estConst.dragCoefficientAir;
    Cw = estConst.windVel;
    Cr = estConst.rudderCoefficient;
    ut = u(1);
    ur = u(2);

    q = [sx;
         sy;
         cos(phi)*(tanh(ut) - Cdh*(sx^2 + sy^2)) - Cda*(sx - Cw*cos(rho))*((sx - Cw*cos(rho))^2 + (sy - Cw*sin(rho))^2)^(1/2);
         sin(phi)*(tanh(ut) - Cdh*(sx^2 + sy^2)) - Cda*(sy - Cw*sin(rho))*((sx - Cw*cos(rho))^2 + (sy - Cw*sin(rho))^2)^(1/2);
         Cr*ur;
         0;
         0];
end


function h = evalh(x, estConst)
    px = x(1);
    py = x(2);
    phi = x(5);
    b = x(7);
    xa = estConst.pos_radioA(1);
    ya = estConst.pos_radioA(2);
    xb = estConst.pos_radioB(1);
    yb = estConst.pos_radioB(2);
    xc = estConst.pos_radioC(1);
    yc = estConst.pos_radioC(2);

    h = [((px - xa)^2 + (py - ya)^2)^(1/2);
         ((px - xb)^2 + (py - yb)^2)^(1/2);
         ((px - xc)^2 + (py - yc)^2)^(1/2);
         b + phi;
         phi];
end



function A = evalA(x, u, estConst)
    sx = x(3);
    sy = x(4); 
    phi = x(5);
    rho = x(6);
    Cdh = estConst.dragCoefficientHydr;
    Cda = estConst.dragCoefficientAir;
    Cw = estConst.windVel;
    ut = u(1);
    
    A = [0, 0,                                                                                                                                                                                       1,                                                                                                                                                                                       0,                                        0,                                                                                                                                                                                 0, 0;
         0, 0,                                                                                                                                                                                       0,                                                                                                                                                                                       1,                                        0,                                                                                                                                                                                 0, 0;
         0, 0, - Cda*((sx - Cw*cos(rho))^2 + (sy - Cw*sin(rho))^2)^(1/2) - 2*Cdh*sx*cos(phi) - (Cda*(sx - Cw*cos(rho))*(2*sx - 2*Cw*cos(rho)))/(2*((sx - Cw*cos(rho))^2 + (sy - Cw*sin(rho))^2)^(1/2)),                                                           - 2*Cdh*sy*cos(phi) - (Cda*(sx - Cw*cos(rho))*(2*sy - 2*Cw*sin(rho)))/(2*((sx - Cw*cos(rho))^2 + (sy - Cw*sin(rho))^2)^(1/2)), -sin(phi)*(tanh(ut) - Cdh*(sx^2 + sy^2)), (Cda*Cw*(sx - Cw*cos(rho))*(sy*cos(rho) - sx*sin(rho)))/((sx - Cw*cos(rho))^2 + (sy - Cw*sin(rho))^2)^(1/2) - Cda*Cw*sin(rho)*((sx - Cw*cos(rho))^2 + (sy - Cw*sin(rho))^2)^(1/2), 0;
         0, 0,                                                           - 2*Cdh*sx*sin(phi) - (Cda*(sy - Cw*sin(rho))*(2*sx - 2*Cw*cos(rho)))/(2*((sx - Cw*cos(rho))^2 + (sy - Cw*sin(rho))^2)^(1/2)), - Cda*((sx - Cw*cos(rho))^2 + (sy - Cw*sin(rho))^2)^(1/2) - 2*Cdh*sy*sin(phi) - (Cda*(sy - Cw*sin(rho))*(2*sy - 2*Cw*sin(rho)))/(2*((sx - Cw*cos(rho))^2 + (sy - Cw*sin(rho))^2)^(1/2)),  cos(phi)*(tanh(ut) - Cdh*(sx^2 + sy^2)), Cda*Cw*cos(rho)*((sx - Cw*cos(rho))^2 + (sy - Cw*sin(rho))^2)^(1/2) + (Cda*Cw*(sy - Cw*sin(rho))*(sy*cos(rho) - sx*sin(rho)))/((sx - Cw*cos(rho))^2 + (sy - Cw*sin(rho))^2)^(1/2), 0;
         0, 0,                                                                                                                                                                                       0,                                                                                                                                                                                       0,                                        0,                                                                                                                                                                                 0, 0;
         0, 0,                                                                                                                                                                                       0,                                                                                                                                                                                       0,                                        0,                                                                                                                                                                                 0, 0;
         0, 0,                                                                                                                                                                                       0,                                                                                                                                                                                       0,                                        0,                                                                                                                                                                                 0, 0];

end



function L = evalL(x, u, estConst)
    sx = x(3);
    sy = x(4); 
    phi = x(5);
    Cdh = estConst.dragCoefficientHydr;
    Cr = estConst.rudderCoefficient;
    ur = u(2);

    L = [                      0,     0, 0, 0;
                               0,     0, 0, 0;
     -Cdh*cos(phi)*(sx^2 + sy^2),     0, 0, 0;
     -Cdh*sin(phi)*(sx^2 + sy^2),     0, 0, 0;
                               0, Cr*ur, 0, 0;
                               0,     0, 1, 0
                               0,     0, 0, 1];
end


function Pdot = evalPdot(P, x, u, estConst)
    A = evalA(x, u, estConst);
    L = evalL(x, u, estConst);
    Qc = diag([estConst.DragNoise, estConst.RudderNoise,...
               estConst.WindAngleNoise, estConst.GyroDriftNoise]);
           
    Pdot = A * P + P * A' + L * Qc * L';
end


function qAndPdot = evalqAndPdot(xAndP, u, estConst, n_states)
    % stack evaluations of q and Pdot for integration in ode45 solver (to compute x_prior and P_prior)
    x = xAndP(1:n_states);
    P = reshape(xAndP(n_states+1:end), n_states, n_states);
    
    q = evalq(x, u, estConst);
    Pdot = evalPdot(P, x, u, estConst);
    
    qAndPdot = [q; reshape(Pdot, [], 1)];
end


function H = evalH(x, estConst)
    px = x(1);
    py = x(2);
    xa = estConst.pos_radioA(1);
    ya = estConst.pos_radioA(2);
    xb = estConst.pos_radioB(1);
    yb = estConst.pos_radioB(2);
    xc = estConst.pos_radioC(1);
    yc = estConst.pos_radioC(2);

    H = [(px - xa)/((px - xa)^2 + (py - ya)^2)^(1/2), (py - ya)/((px - xa)^2 + (py - ya)^2)^(1/2), 0, 0, 0, 0, 0;
         (px - xb)/((px - xb)^2 + (py - yb)^2)^(1/2), (py - yb)/((px - xb)^2 + (py - yb)^2)^(1/2), 0, 0, 0, 0, 0;
         (px - xc)/((px - xc)^2 + (py - yc)^2)^(1/2), (py - yc)/((px - xc)^2 + (py - yc)^2)^(1/2), 0, 0, 0, 0, 0;
                                                   0,                                           0, 0, 0, 1, 0, 1;
                                                   0,                                           0, 0, 0, 1, 0, 0];

end

function M = evalM()
    M = eye(5);
end