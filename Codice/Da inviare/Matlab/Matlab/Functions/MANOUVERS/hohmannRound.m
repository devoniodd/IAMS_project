function [dV,tof,orbT] = hohmannRound(orb1,orb2)
% hohmannRound - Hohmann maneuver between round orbits
%
% PROTOTYPE:
% [dV,tof,orbT] = hohmannRound(orb1,orb2)
%
% DESCRIPTION:
% This function computes the transfer orbit between two circular orbits
% given their orbital parameters using Hohmann's maneuver and outputs the 
% delta velocities required, the time of flight and the orbital parameters 
% of the transfer orbit.
%
% INPUT:
% initial orbit     [1x7]           Initial orbital parameters                          [N/D]
% final orbit       [1x7]           Final orbital parameters                            [N/D]
%
% OUTPUT:
% delta V           [2x1]           Vector containg the two delta V required            [km/s]
% time of flight    [1x1]           Time to tranfer from orbit 1 to orbit 2             [s]
% transfer orbit    [1x7]           Elliptical transfer orbit parameters                [N/D]
%                                                       
%% UTILS
if ismac
    load("../Data/utils.mat",'mu');
else
    load("..\Data\utils.mat",'mu');
end

%% VARIABLES CHECK
if orb1(1,3) ~= orb2(1,3)
    error('The orbits must belong to the same plane')
end

%% INPUT EXTRACTION

% Radius
r1 = orb1(1,1);
r2 = orb2(1,1);

%Velocities
V1 = sqrt(mu/r1);
V2 = sqrt(mu/r2);

%% TRANSFER ORBIT CHARACTERIZATION

% Vector creation
orbT = zeros(1,7);

% Orbital parameters
rTa = r2;
rTp = r1;
orbT(1,1) = (rTa+rTp)/2;
orbT(1,2) = (rTa-rTp)/(rTa+rTp);
orbT(1,3) = orb1(1,3);
orbT(1,4) = orb1(1,4);
orbT(1,5) = orb1(1,5) + orb1(1,6);
orbT(1,6) = 0;
orbT(1,7) = pi;
pT = orbT(1,1)*(1-orbT(1,2)^2);

% Velocities
VTp = sqrt(mu/pT)*(1+orbT(1,2));
VTa = sqrt(mu/pT)*(1-orbT(1,2));


%% TRANSFER

% Delta Vs
dV1 = VTp-V1;
dV2 = V2-VTa;

dV = [dV1; dV2];

% Tranfer Time 
tof = pi * sqrt((orbT(1,1)^3)/mu);


end