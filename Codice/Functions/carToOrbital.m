function [a,e,i,Omega,omega,theta] = carToOrbital(r,v)
% rvToOrbit - Conversion from Cartesian coordinates to Keplerian elements
%
% PROTOTYPE:
% [a,e,i,Omega,omega,theta] = rvToOrbital(r,v)
%
% DESCRIPTION:
% Conversion from Cartesian coordinates to Keplerian elements. Angles in
% radians.
%
% INPUT:
% r                 [3x1]           Position vector                 [km]
% v                 [3x1]           Velocity vector                 [km/s]
%
% OUTPUT:
% a                 [1x1]           Semi-major axis                 [km]
% e                 [1x1]           Eccentricity                    [-]
% i                 [1x1]           Inclination                     [rad]
% Omega             [1x1]           RAAN                            [rad]
% omega             [1x1]           Pericentre anomaly              [rad]
% theta             [1x1]           True anomaly                    [rad]


%% VALUE CHECK
if nargin < 2
    error("Error: not enough variables")
end

if length(r) ~= 3 || length(v) ~= 3
    error("Please provide a 1x3 vector for both radius and velocity");
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

%% SEMI - MAJOR AXIS
a = 1/(2/rNorm - vNorm^2/mu);

%% ANGULAR MOMENTUM
h = cross(r,v);
hNorm = norm(h);

%% ECCENTRICITY
e = 1/mu * (cross(v,h) - mu*r/rNorm);
eNorm = norm(e);

%% INCLINATION
i = acos(dot(h,kDir)/hNorm);

%% NODAL PLANE DIRECTION
N = cross(kDir,h);
NNorm = norm(N);

%% RIGHT ASCENSION OF THE ASCENDING NODE - LONGITUDE OF ASCENDING NODE
if dot(N,jDir) >= 0
    Omega = acos( dot(N,iDir) / NNorm);
else 
    Omega = 2*pi - acos( dot(N,iDir) / NNorm);
end

%% ARGUMENT OF PERIAPSIS
if dot(e,kDir) >= 0 
    omega = acos(dot(N,e) / (NNorm * eNorm));
else 
    omega = 2*pi - acos(dot(N,e) / (NNorm * eNorm));
end

%% MEAN ANOMALY

if dot(v,r) >= 0
    theta = acos(dot(e,r) / (eNorm * rNorm));
else 
    theta = 2*pi - acos(dot(e,r) / (eNorm * rNorm));
end


end