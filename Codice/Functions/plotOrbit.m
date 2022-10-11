function [X,Y,Z,V] = plotOrbit(kepEl, stepTh, degrees, deltaTh, startTh)
% plotOrbit.m - Plot the arc length deltaTh of the orbit described by kepEl.
%
% \\\\\\\\\\\\\\\\\\\\\\\\\\
%       TO BE EDITED
% \\\\\\\\\\\\\\\\\\\\\\\\\\
%
% PROTOTYPE:
% plotOrbit(kepEl, mu, deltaTh, stepTh)
%
% DESCRIPTION:
% Plot the arc length of the orbit described by a set of orbital
% elements for a specific arc length.
%
% INPUT:
% kepEl [1x6] orbital elements [km,rad]
% mu [1x1] gravitational parameter [km^3/s^2]
% deltaTh [1x1] arc length [rad]
% stepTh [1x1] arc length step [rad]
%
% OUTPUT:
% X [1xn] X position [km]
% Y [1xn] Y position [km]
% Z [1xn] Z position [km]
% V [1xn] V norm     [km/s]


%% CHECK OF VARIABLES
if length(kepEl) < 1
    error("Error: please insert valid parameters");
end

%% EXTRACTING VARIABLES

a = kepEl(1);
e = kepEl(2);
i = kepEl(3);
Omega = kepEl(4);
omega = kepEl(5);
theta = kepEl(6);

%% DEG TO RAD
if nargin == 3 && degrees ~= 1
    if degrees ~= 0
        error(sprintf("Please select a valid option: \n1 - Input in degrees \n0 - Input in radians"));
    end
end

if degrees == 1
    Omega = deg2rad(Omega);
    omega = deg2rad(omega);
    theta = deg2rad(theta);
    i = deg2rad(i);
end

%% VALUE CHECK
% To be added - Deg2Rad conversion!!
if nargin == 5 && startTh < theta
    error("Error: plotted orbit doesn't starts after theta");
end

%% INITIAL CONDITIONS

% Orbit arc starts from theta
if nargin < 5
    startTh = theta;
end

% Orbit arc is 360 deg
if nargin < 4
    deltaTh = 2*pi;
end

% stepTh
if nargin < 2
    stepTh = pi/180;
end

%% PLOT ORBIT
thetaV = (startTh : stepTh : startTh+deltaTh); % Declaration of vector theta
stepTothetaV = (startTh : stepTh : theta);
stepTotheta = length(stepTothetaV);
Nstep = length(thetaV);

X = zeros(1,Nstep);
Y = zeros(1,Nstep);
Z = zeros(1,Nstep);
V = zeros(1,Nstep);

%% ORBIT TO RV
for j = 1 : Nstep
    [r,v] = orbitalToCar(a,e,i,Omega,omega,thetaV(j));
    X(1,j) = r(1);
    Y(1,j) = r(2);
    Z(1,j) = r(3);
    V(j) = norm(v);
end

end


