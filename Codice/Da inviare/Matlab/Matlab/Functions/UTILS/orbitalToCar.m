function [r,v] = orbitalToCar(a,e,i,Omega,omega,theta)
% Conversion from Keplerian elements to Cartesian coordinates

% INPUT:
% a         [1x1]   Semi-major axis [km]
% e         [1x1]   Eccentricity [-]
% i         [1x1]   Inclination, if degrees = 0 [rad] else if degrees = 1 [deg]
% Omega     [1x1]   RAAN, if degrees = 0 [rad] else if degrees = 1 [deg]
% omega     [1x1]   Pericentre anomaly, if degrees = 0 [rad] else if degrees = 1 [deg]
% theta     [1x1]   True anomaly, if degrees = 0 [rad] else if degrees = 1 [deg]
%
% OUTPUT:
% r         [3x1]   Position vector [km]
% v         [3x1]   Velocity vector [km/s]

%% VALUE CHECK
if nargin < 6
    error("Please insert a valid amount of variables");
end

if length(e) == 3 
    e = norm(e);
elseif length(e) ~= 1
    error("Please check e value, it is neither a 3 element vector nor a scalar!");
end

%% UTILS IMPORT 
if ismac
    load("../Data/utils.mat",'mu');
else
    load("..\Data\utils.mat",'mu');
end

%% R,V PERIFOCAL SYSTEM
if e >= 0 && e < 1 % Elliptic
    p = a*(1-e^2); % Semilat rett
elseif e == 1 % Parabolic
    p = 2 * a;
else
    p = a*(e^2 - 1);    
end

rNorm = p/(1+e*cos(theta)); % Radius
rPF = rNorm * [cos(theta); sin(theta); 0];
vPF = sqrt(mu/p) * [-sin(theta); (e+cos(theta)); 0];
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

