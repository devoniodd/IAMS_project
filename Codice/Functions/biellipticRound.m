function [dV,tof,orbT1,orbT2] = biellipticRound(orb1,orb2,rMax)
% biellipticRound - Bielliptic maneuver between round orbits
%
% PROTOTYPE:
% [dV,tof,orbT1,orbT2] = biellipticRound(orb1,orb2,rMax)
%
% DESCRIPTION:
% This function computes the transfer orbit between two circular orbits
% given their orbital parameters and the maximum radius using the bielliptic
% maneuver and outputs the delta velocities required, the time of flight and 
% the orbital parameters of the two transfer orbits.
%
% INPUT:
% initial orbit         [1x7]           Initial orbital parameters                          [N/D]
% final orbit           [1x7]           Final orbital parameters                            [N/D]
% maximum radius        [1x1]           Maximum radius used in the transfer orbits          [Km]
%
% OUTPUT:
% delta V               [3x1]           Vector containg the three delta V required          [km/s]
% time of flight        [2x1]           Vector containg the time of flight for each orbit   [s]
% First transfer orbit  [1x7]           First transfer orbit parameters                     [N/D]
% Second transfer orbit [1x7]           second transfer orbit parameters                    [N/D]
%   
%% UTILS
if ismac
    load("../Data/utils.mat",'mu');
else
    load("..\Data\utils.mat",'mu');
end

%% VARIABLES CHECK
if rMax < orb2(1,1)
    error("The maximum radius can't be smaller tha the radius of the second orbit")
end

if orb1(1,3) ~= orb2(1,3)
    error('The orbits must belong to the same plane')
end

%% INITIAL AND FINAL ORBITS

% Radius extraction
ri = orb1(1,1);
rf = orb2(1,1);

% Velocities
Vi = sqrt(mu/ri);
Vf = sqrt(mu/rf);

%% TRANSFER ORBIT 1 

% Vector creation
orbT1 = zeros(1,7);

% Orbital parameters

orbT1(1,1) = (ri+rMax)/2;
orbT1(1,2) = (rMax-ri)/(rMax+ri);
orbT1(1,3) = orb1(1,3);
orbT1(1,4) = orb1(1,4);
orbT1(1,5) = orb1(1,5) + orb1(1,6);
orbT1(1,6) = 0;
orbT1(1,7) = pi;
pT1 = orbT1(1,1)*(1-orbT1(1,2)^2);

% Velocities
VT1p = sqrt(mu/pT1)*(1+orbT1(1,2));
VT1a = sqrt(mu/pT1)*(1-orbT1(1,2));

%% TRANSFER ORBIT 2 

% Vector creation
orbT2 = zeros(1,7);

% Orbital parameters

orbT2(1,1) = (rf+rMax)/2;
orbT2(1,2) = (rMax-rf)/(rMax+rf);
orbT2(1,3) = orb2(1,3);
orbT2(1,4) = orb2(1,4);
orbT2(1,5) = orb2(1,5) + orb2(1,6);
orbT2(1,6) = pi;
orbT2(1,7) = 0;
pT2 = orbT2(1,1)*(1-orbT2(1,2)^2);

% Velocities
VT2p = sqrt(mu/pT2)*(1+orbT2(1,2));
VT2a = sqrt(mu/pT2)*(1-orbT2(1,2));

%% OUTPUTS

% Velocities 
dV1 = VT1p-Vi;
dV2 = VT2a-VT1a;
dV3 = Vf-VT2p;

dV = [dV1; dV2; dV3];

% Time of flight
t1 = pi*sqrt(orbT1(1,1)^3/mu);
t2 = pi*sqrt(orbT2(1,1)^3/mu);

tof = [t1; t2];



