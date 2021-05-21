% Define state
syms px py sx sy phi rho b
% Define input
syms ut ur
% Define measurements
syms za zb zc zg zn
% Define process noise
syms vd vr vrho vb
% Define measurement noise
syms wa wb wc wg wn
% Define drag constants
syms Cdh Cda Cw Cr
% Define station positions
syms xa xb xc ya yb yc

% Define vectors for input, state noises
x = [px;
     py;
     sx;
     sy;
     phi;
     rho;
     b];

u = [ut;
     ur];
  
v = [vd;
     vr;
     vrho;
     vb];

z = [za;
     zb;
     zc;
     zg;
     zn];

w = [wa;
     wb;
     wc;
     wg;
     wn];
    
% Define state equation
q = [sx;
     sy;
     cos(phi)*(tanh(ut) - Cdh*(sx^2+sy^2)*(1+vd)) - Cda*(sx - Cw*cos(rho))*sqrt((sx - Cw*cos(rho))^2 + (sy - Cw*sin(rho))^2);
     sin(phi)*(tanh(ut) - Cdh*(sx^2+sy^2)*(1+vd)) - Cda*(sy - Cw*sin(rho))*sqrt((sx - Cw*cos(rho))^2 + (sy - Cw*sin(rho))^2);
     Cr*ur*(1+vr);
     vrho;
     vb];
     
% Define measurement equations
h = [sqrt((px - xa)^2 + (py - ya)^2) + wa;
     sqrt((px - xb)^2 + (py - yb)^2) + wb;
     sqrt((px - xc)^2 + (py - yc)^2) + wc;
     phi + b + wg;
     phi + wn];

% Set v equal to sero to evaluate the jacobians at v=0
vd = 0;
vr = 0;
vrho = 0;
vb = 0;

% Q with v = 0
qv0 = simplify(eval(q));

% Calculate Prior update Jacobians
A = simplify(eval(jacobian(q, x)));
L = simplify(eval(jacobian(q, v)));

% Set w equal to sero to evaluate the jacobians at w=0
wa = 0;
wb = 0;
wc = 0;
wg = 0;
wn = 0;

% h with w=0
hw0 = simplify(eval(h));

% Calculate Posterior update Jacobians
H = simplify(eval(jacobian(h, x)));
M = eval(jacobian(h, w));
     
     
     