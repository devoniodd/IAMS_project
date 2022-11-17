function [r,v] = orbitalToCar(orbit)
% orbitalToCar - Conversion from Keplerian elements to Cartesian coordinates
%
% PROTOTYPE:
% [r,v] = orbitalToCar(orbit)
%
% DESCRIPTION:
% Conversion from Keplerian elements to Cartesian coordinates. 
% Angles must be given in radians.
% The input is a standardized vector which contains all of the orbital parameters.
%
% INPUT:
% orbit     [1x7]   Orbital parameters vector   [N/D]
%
% OUTPUT:
% r         [3x1]   Position vector             [km]
% v         [3x1]   Velocity vector             [km/s]

%% UTILS IMPORT 
if ismac
    load("../Data/utils.mat",'mu');
else
    load("..\Data\utils.mat",'mu');
end

%% VECTOR EXTRACTION
a = orbit(1,1);
e = orbit(1,2);
i = orbit(1,3);
Omega = orbit(1,4);
omega = orbit(1,5);
theta = orbit(1,6);

%% RADIUS CALCULATION

% Semi-latus rectum
p = a*(1-e^2);

% Radius norm
rNorm = p/(1+e*cos(theta));

%% R,V IN PERIFOCAL R.S. 
rPF = rNorm * [cos(theta); sin(theta); 0];
vPF = sqrt(mu/p) * [-sin(theta); (e+cos(theta)); 0];

%% FROM PERIFOCAL R.S. TO GEOCENTRIC EQUATORIAL R.S.

% Rotation of angle Omega along earth's K axis
R1 = [cos(Omega),sin(Omega),0; -sin(Omega),cos(Omega),0; 0,0,1];

% Rotation of angle i along our orbit's I axis
R2 = [1,0,0; 0,cos(i),sin(i); 0,-sin(i),cos(i)];

% Rotation of angle omega along our orbit's K axis
R3 = [cos(omega),sin(omega),0; -sin(omega),cos(omega),0; 0,0,1];

% Final rotation matrix
T = R1'*R2'*R3';

% From Perifocal R.S. to Geocentric Equatorial R.S.
r = T * rPF;
v = T * vPF;

end

