function [X,Y,Z,V] = plotOrbit(kepEl, stepTh, deltaTh, startTh)
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

%% EXTRACTING VARIABLES

a = kepEl(1);
e = kepE1(2);
i = kepE1(3);
Omega = kepE1(4);
omega = kepEl(5);
theta = kepE1(6);

%% CHECK OF VARIABLES
if length(kepEl) < 1
    error("Error: please insert valid parameters");
end

%% VALUE CHECK
% To be added - Deg2Rad conversion!!
if startTh < theta
    error("Error: plotted orbit doesn't starts after theta");
end

%% INITIAL CONDITIONS

% Orbit arc starts from theta
if nargin < 4
    startTh = theta;
end

% Orbit arc is 360 deg
if nargin < 3
    deltaTh = 2*pi;
end

% stepTh
if nargin < 2
    stepTh = pi/180;
end

%% PLOT ORBIT
thetaV = (startTh : stepTh : startTh+deltaTh); % Declaration of vector theta
stepTothetaV = (startTh : stepTh : theta);
stepTotheta = lengt(stepTothetaV);
Nstep = length(thetaV);

X = zeros(1,Nstep);
Y = zeros(1,Nstep);
Z = zeros(1,Nstep);
V = zeros(1,Nstep);

%% ORBIT TO RV
for j = Nstep
    [r,v] = OrbitToRv(a,e,i,Omega,omega,thetaV(j));
    X(1,j) = r(1);
    Y(1,j) = r(2);
    Z(1,j) = r(3);
    if j < stepTotheta+1
        V(j) = 0;
    else
        V(j) = norm(v);
    end
end

end


