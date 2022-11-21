function [dV,phiV] =  apseLineRotation(eta,rai,rpi,raf,rpf)

%% UTILS
if ismac
    load("../Data/utils.mat",'mu');
else
    load("..\Data\utils.mat",'mu');
end

%% ORBITS DEFINITION
ei = (rai-rpi)/(rai+rpi);
ef = (raf-rpf)/(raf+rpf);

hi = sqrt(rpi * (1+ei) * mu);
hf = sqrt(rpf * (1+ef) * mu);

%% INTERSECTION POINTS

% Problem definition
a = ei*hf^2 - ef*hi^2 * cos(eta);
b = -ef*hi^2 * sin(eta);
c = hi^2 - hf^2;
phi = atan(b/a);

% Intersection point anomalies
thetaI = phi + acos(c/a * cos(phi));
thetaJ = phi - acos(c/a * cos(phi));

if thetaI < 0
    thetaI = 2*pi-thetaI;
end
if thetaJ < 0
    thetaJ = 2*pi-thetaJ;
end

%% MANOUVER POINT

% Anomaly
if thetaI>thetaJ
    theta = thetaJ;
else 
    theta = thetaI;
end

% Radius
r = (hi^2)/mu * 1/(1+ei*cos(theta));

% Velocities
Vti = hi/r;
Vri = mu/hi * ei * sin(theta);
Vi = sqrt(Vti^2 + Vri^2);

Vtf = hf/r;
Vrf = mu/hf * ef * sin(theta-eta);
Vf = sqrt(Vtf^2 + Vrf^2);

% Flight path angle
gammai = atan(Vri/Vti);
gammaf = atan(Vrf/Vtf);

%% OUTPUTS

dV = sqrt(Vi^2 + Vf^2 - 2*Vi*Vf*cos(gammaf-gammai));
phiV = atan((Vrf-Vri)/(Vtf-Vti));
if phiV < 0
    phiV = phiV + pi;
end






