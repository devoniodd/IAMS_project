function [r,v] = orbitalToCar(a,e,i,Omega,omega,theta,degrees)
% Conversion from Keplerian elements to Cartesian coordinates

% INPUT:
% a         [1x1]   Semi-major axis [km]
% e         [1x1]   Eccentricity [-]
% i         [1x1]   Inclination, if degrees = 0 [rad] else if degrees = 1 [deg]
% Omega     [1x1]   RAAN, if degrees = 0 [rad] else if degrees = 1 [deg]
% omega     [1x1]   Pericentre anomaly, if degrees = 0 [rad] else if degrees = 1 [deg]
% theta     [1x1]   True anomaly, if degrees = 0 [rad] else if degrees = 1 [deg]
% mu        [1x1]   Gravitational parameter [km^3/s^2]
%
% OUTPUT:
% r         [3x1]   Position vector [km]
% v         [3x1]   Velocity vector [km/s]

%% VALUE CHECK
if nargin < 6
    error("Please insert a valid amount of variables");
end

if nargin == 7 && degrees ~= 1
    if degrees ~= 0
        error(sprintf("Please select a valid option: \n1 - Input in degrees \n0 - Input in radians"));
    end
%    disp("Input is in degrees");
else
%    disp("Input is in radians");
    degrees = 0;
end

%% DEG TO RAD
if degrees == 1
    Omega = deg2rad(Omega);
    omega = deg2rad(omega);
    theta = deg2rad(theta);
    i = deg2rad(i);
end

%% UTILS IMPORT 
if ismac
    load("../Data/utils.mat",'mu');
else
    load("..\Data\utils.mat",'mu');
end

%% VECTORS NORMALIZATION
eNorm = norm(e);

%% CALCOLI ...sistema titolo

% Semilat rett
p = a*(1-eNorm^2);

% Radius
rNorm = p/(1+eNorm*cos(theta));

%% R,V PERIFOCAL SYSTEM
rPF = rNorm * [cos(theta); sin(theta); 0];
vPF = sqrt(mu/p) * [-sin(theta); (eNorm+cos(theta)); 0];

%% FROM PERIFOCAL TO GEOCENTRIC

% Rotation of angle Omega along earth's K axis
R1 = [cos(Omega),sin(Omega),0; -sin(Omega),cos(Omega),0; 0,0,1];

% Rotation of angle i along our orbit's I axis
R2 = [1,0,0; 0,cos(i),sin(i); 0,-sin(i),cos(i)];

% Rotation of angle omega along our orbit's K axis
R3 = [cos(omega),sin(omega),0; -sin(omega),cos(omega),0; 0,0,1];

% Rotation
T = R1'*R2'*R3';

% From Perifocal to Geocentric Equatorial
r = T * rPF;
v = T * vPF;

end

