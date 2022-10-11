function [dV,tT,tOrbit] = hohmannRound(r1,r2,dVsum)
% This function computes the transfer orbit between two circular orbits
% INPUT:
% r1 - The radius of the starting circular orbit
% r2 - The radius of the target circular orbit
% OUTPUT:
% tOrbit - Structure containing all parameters regarding the elliptical
% transfer orbit
% dV - Sum of the delta Vs for the two impulses [if dVsum is 1] or a vector
% containing two values representing the two impulses [id dVsum is 0]
% tT - Time to tranfer from orbit 1 to orbit 2

%% UTILS
if ismac
    load("../Data/utils.mat",'mu');
else
    load("..\Data\utils.mat",'mu');
end

%% TRANSFER ORBIT CHARACTERIZATION
% Orbital parameters
rTa = r2;
rTp = r1;
eT = (rTa-rTp)/(rTa+rTp);
aT = (rTa+rTp)/2;
pT = aT*(1-eT^2);

% Velocities
VTp = sqrt(mu/pT)*(1+eT);
VTa = sqrt(mu/pT)*(1-eT);

V1 = sqrt(mu/r1);
V2 = sqrt(mu/r2);

%% TRANSFER

% Delta Vs
dV1 = VTp-V1;
dV2 = V2-VTa;

% Tranfer Time 
tT = pi * sqrt((aT^3)/mu);

%% OUTPUTS

if dVsum
    dV = dV1 + dV2;
else
    dV = [dV1, dV2];
end

tOrbit.e = eT;
tOrbit.a = aT;
tOrbit.p = pT;
tOrbit.ra = rTa;
tOrbit.rp = rTp;
tOrbit.Vp = VTp;
tOrbit.Va = VTa;

end