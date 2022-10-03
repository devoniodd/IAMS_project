function [a,e,i,Omega,omega,theta] = rvToOrbital(r,v,mu,degrees)
% Description - TBA

%% VALUE CHECK
if nargin <= 2
    error("Error: not enough variables")
end

if nargin == 4 && degrees ~= "Yes"
    if degrees ~= "No"
        error("Please select a valid option:\nYes - Output in degrees\nNo - output in radians");
    end
    disp("Answer is in radians");
    deg = 0;
else
    disp("Answer is in degrees");
    deg = 1;
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
a = (2/rNorm - vNorm^2/mu).^(-1);

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
    Omega = acos(nNorm*iDir./nNorm);
else 
    Omega = 2*pi - acos(nNorm*iDir./nNorm);
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
if nargin == 3 || deg == 1
    Omega = 0.5 * Omega/pi * 360;
    omega = 0.5 * omega/pi * 360;
    theta = 0.5 * theta/pi * 360;
end

end