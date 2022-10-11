function [a,e,i,Omega,omega,theta] = carToOrbital(r,v,degrees)
% rvToOrbit - Conversion from Cartesian coordinates to Keplerian elements
%
% PROTOTYPE:
% [a,e,i,Omega,omega,theta] = rvToOrbital(r,v,mu,degrees)
%
% DESCRIPTION:
% Conversion from Cartesian coordinates to Keplerian elements. Angles in
% radians.
%
% INPUT:
% r                 [3x1]           Position vector                 [km]
% v                 [3x1]           Velocity vector               [km/s]
% mu                [1x1]           Gravitational parameter   [km^3/s^2]
%
% OUTPUT:
% a                 [1x1]           Semi-major axis                 [km]
% e                 [1x1]           Eccentricity                     [-]
% i                 [1x1]           Inclination                [rad/deg]*
% Omega             [1x1]           RAAN                       [rad/deg]*
% omega             [1x1]           Pericentre anomaly         [rad/deg]*
% theta             [1x1]           True anomaly               [rad/deg]*
%
% * Default value is in radians

%% VALUE CHECK
if nargin <= 2
    error("Error: not enough variables")
end

if nargin == 4 && degrees ~= 1
    if degrees ~= 0
        error("Please select a valid option:    1 - Output in degrees      0 - output in radians");
    end
    disp("Answer is in radians");
else
    disp("Answer is in degrees");
end

%% UTILS IMPORT 
if ismac
    load("../Data/utils.mat",'mu','iDir','jDir','kDir');
else
    load("..\Data\utils.mat",'mu','iDir','jDir','kDir');
end

%% VECTOR NORMALIZATION
rNorm = norm(r);
vNorm = norm(v);
% vRadial = r*v/rNorm; % Non necessari

%% SEMI - MAJOR AXIS
a = 1/(2/rNorm - vNorm^2/mu);

%% ANGULAR MOMENTUM
h = cross(r,v);
hNorm = norm(h);

%% ECCENTRICITY
e = 1/mu * (cross(v,h) - mu*r./rNorm);
eNorm = norm(e);

%% INCLINATION
i = acos(dot(h,kDir)/hNorm);

%% NODAL PLANE
n = cross(kDir,h) / norm (cross(kDir,h));
nNorm = norm(n);

%% RIGHT ASCENSION OF THE ASCENDING NODE - LONGITUDE OF ASCENDING NODE
if nNorm*jDir >= 0 % Condizioni da sistemare
    Omega = acos( dot(n,iDir) / nNorm);
else 
    Omega = 2*pi - acos( dot(n,iDir) / nNorm);
end

%% ARGUMENT OF PERIAPSIS
if eNorm*nNorm >= 0 % Ste condizioni sono da sistemare
    omega = acos(dot(e,n) / eNorm);
else 
    omega = 2*pi - acos(dot(e,n) / eNorm);
end

%% MEAN ANOMALY

if dot(v,r) >= 0  % Anche queste porcozzi
    theta = acos(dot(e/eNorm,r/rNorm));
else 
    theta = 2*pi - acos(dot(e/eNorm,r/rNorm));
end

%% RADINATS TO DEGREES
if nargin == 3 && degrees == 1
    Omega = rad2deg(Omega);
    omega = rad2deg(omega);
    theta = rad2deg(theta);
    i = rad2deg(i);
end

end