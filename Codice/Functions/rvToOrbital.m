function [a,e,i,Omega,omega,theta] = rvToOrbital(r,v,mu,degrees)
% Description - TBA

%% VALUE CHECK
if nargin <= 2
    error("Error: not enough variables")
end

if nargin == 4 && degrees ~= 1
    if degrees ~= 0
        error("Please select a valid option:\n1 - Output in degrees\n0 - output in radians");
    end
    disp("Answer is in radians");
else
    disp("Answer is in degrees");
end

%% DIRECTIONAL VERSORS
iDir = [1; 0; 0];
jDir = [0; 1; 0];
kDir = [0; 0; 1];

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
i = acos(hNorm*kDir/hNorm);

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
if rNorm*eNorm >= 0  % Anche queste porcozzi
    theta = acos(dot(r,n) / (rNorm*eNorm));
else 
    theta = 2*pi - acos(dot(r,n) / (rNorm*eNorm));
end

%% RADINATS TO DEGREES
if nargin == 3 || degrees == 1
    Omega = rad2deg(Omega);
    omega = rad2deg(omega);
    theta = rad2deg(theta);
end

end