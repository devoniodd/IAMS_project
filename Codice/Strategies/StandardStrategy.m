clear 
close all
clc

%% STANDARD STRATEGY

% This matlab script shows and calculates the standard strategy to transfer
% between two given orbits. The strategy is defined as follows:
% - Change of plane;
% - Change of Periapsis argument;
% - Bitangent transfer starting from the periapsis and arriving at the apoapsis
% 
% At the end you will find all the orbital parameters of the orbits used, a
% vector containing all the velocities differences, the total velocity
% difference required, a vector containing all the time required for each
% manuever and the total time of flight for the entire transfer

%% UTILS
if ismac
    load("../Data/utils.mat",'mu');
else
    load("..\Data\utils.mat",'mu');
end

%% INITIAL AND FINAL ORBITS

% Initial Radius 
r = [-8048.2861; -4171.3048; 2895.9296];

% Initial Velocity
V = [1.7540; -5.9910; -1.9520]; 

% Initial Orbit 
Orbit1 = zeros(1,7);
[a,e,i,O,o,thi] = carToOrbital(r,V);
Orbit1(1,1) = a;
Orbit1(1,2) = e;
Orbit1(1,3) = i;
Orbit1(1,4) = O;
Orbit1(1,5) = o;
Orbit1(1,6) = thi;
Orbit1(1,7) = nan;

CurrentOrbit = Orbit1;

% Final Orbit
OrbitFinal = zeros(1,7);
OrbitFinal(1,1) = 16720;
OrbitFinal(1,2) = 0.2502;
OrbitFinal(1,3) = 1.1190;
OrbitFinal(1,4) = 0.6245;
OrbitFinal(1,5) = 3.1350;
OrbitFinal(1,6) = nan;
OrbitFinal(1,7) = 3.1000;


%% PLANE CHANGE

[dv1,t1,Orbit1,Orbit2] = planeChange(CurrentOrbit,OrbitFinal(1,3),OrbitFinal(1,4));

dV = [dv1];
t = [t1];
CurrentOrbit = Orbit2;

%% PERIAPSIS CHANGE

[dv2,t2,Orbit2,Orbit3] = changePeriapsisArg(CurrentOrbit,OrbitFinal(1,5));

dV = [dV; dv2];
t = [t; t2];
CurrentOrbit = Orbit3;

%% SHAPE AND DIMENSION CHANGE

% Maneuver Point
Orbit3(1,7) = 0;
CurrentOrbit(1,7) = 0;

% Time to reach maneuver point
t3 = timeOfFlight(CurrentOrbit,CurrentOrbit(1,6),CurrentOrbit(1,7));
t = [t; t3];

% Transfer Orbit
rtp = CurrentOrbit(1,1)*(1-CurrentOrbit(1,2));
rta = OrbitFinal(1,1)*(1+OrbitFinal(1,2));

TransferOrbit = zeros(1,7);
TransferOrbit(1,1) = (rtp+rta)/2;
TransferOrbit(1,2) = (rta-rtp)/(rta+rtp);
TransferOrbit(1,3) = CurrentOrbit(1,3);
TransferOrbit(1,4) = CurrentOrbit(1,4);
TransferOrbit(1,5) = CurrentOrbit(1,5);
TransferOrbit(1,6) = 0;
TransferOrbit(1,7) = pi;

Orbit4 = TransferOrbit;

% Velocities 
pCurrent = CurrentOrbit(1,1)*(1-CurrentOrbit(1,2)^2);
pFinal = OrbitFinal(1,1)*(1-OrbitFinal(1,2)^2);
PTransfer = (2*rtp*rta)/(rtp+rta);

Vpc = sqrt(mu/pCurrent)*(1+CurrentOrbit(1,2));
Vaf = sqrt(mu/pFinal)*(1-OrbitFinal(1,2));
Vpt = sqrt(mu/PTransfer)*(1+TransferOrbit(1,2));
Vat = sqrt(mu/PTransfer)*(1-TransferOrbit(1,2));

dv3 = Vpt - Vpc;
dv4 = Vaf - Vat;
dV = [dV; dv3; dv4];

% Transfer orbit time of flight
t4 = pi*sqrt(TransferOrbit(1,1)^3/mu);
t = [t; t4];

% Time to reach destination 
OrbitFinal(1,6) = pi;
CurrentOrbit = OrbitFinal;

t5 = timeOfFlight(CurrentOrbit,CurrentOrbit(1,6),CurrentOrbit(1,7));
t = [t; t5];

%% FINAL OUTPUTS

TotalT = sum(t);
TotalDV = sum(abs(dV));

%% PLOT

OrbitsToPlot = [Orbit1; Orbit2; Orbit3; Orbit4; OrbitFinal];

drawOrbitapp(OrbitsToPlot,0.01,0)