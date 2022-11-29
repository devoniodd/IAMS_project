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

%% DATA AND UTILS IMPORT 

addpath(genpath("../Data/"))
load("GivenData.mat");
load("utils.mat",'mu');

%% INITIAL AND FINAL ORBITS

% Initial Orbit 
Orbit1 = zeros(1,7);
[a,e,i,O,o,thi] = carToOrbital(r1,V1);
Orbit1(1) = a;
Orbit1(2) = e;
Orbit1(3) = i;
Orbit1(4) = O;
Orbit1(5) = o;
Orbit1(6) = thi;
Orbit1(7) = nan;

CurrentOrbit = Orbit1;

%% PLANE CHANGE

[dv1,t1,Orbit1,Orbit2] = planeChange(CurrentOrbit,FinalOrbit(3),FinalOrbit(4));

dV = [dv1];
t = [t1];
CurrentOrbit = Orbit2;

%% PERIAPSIS CHANGE

[dv2,t2,Orbit2,Orbit3] = changePeriapsisArg(CurrentOrbit,FinalOrbit(5));

dV = [dV; dv2];
t = [t; t2];
CurrentOrbit = Orbit3;

%% SHAPE AND DIMENSION CHANGE

% Maneuver Point
Orbit3(7) = 0;
CurrentOrbit(7) = 0;

% Time to reach maneuver point
t3 = timeOfFlight(CurrentOrbit,CurrentOrbit(6),CurrentOrbit(7));
t = [t; t3];

% Transfer Orbit
rtp = CurrentOrbit(1)*(1-CurrentOrbit(2));
rta = FinalOrbit(1)*(1+FinalOrbit(2));

TransferOrbit = zeros(1,7);
TransferOrbit(1) = (rtp+rta)/2;
TransferOrbit(2) = (rta-rtp)/(rta+rtp);
TransferOrbit(3) = CurrentOrbit(3);
TransferOrbit(4) = CurrentOrbit(4);
TransferOrbit(5) = CurrentOrbit(5);
TransferOrbit(6) = 0;
TransferOrbit(7) = pi;

Orbit4 = TransferOrbit;

% Velocities 
pCurrent = CurrentOrbit(1)*(1-CurrentOrbit(2)^2);
pFinal = FinalOrbit(1)*(1-FinalOrbit(2)^2);
PTransfer = (2*rtp*rta)/(rtp+rta);

Vpc = sqrt(mu/pCurrent)*(1+CurrentOrbit(2));
Vaf = sqrt(mu/pFinal)*(1-FinalOrbit(2));
Vpt = sqrt(mu/PTransfer)*(1+TransferOrbit(2));
Vat = sqrt(mu/PTransfer)*(1-TransferOrbit(2));

dv3 = Vpt - Vpc;
dv4 = Vaf - Vat;
dV = [dV; dv3; dv4];

% Transfer orbit time of flight
t4 = pi*sqrt(TransferOrbit(1)^3/mu);
t = [t; t4];

% Time to reach destination 
FinalOrbit(6) = pi;
CurrentOrbit = FinalOrbit;

t5 = timeOfFlight(CurrentOrbit,CurrentOrbit(6),CurrentOrbit(7));
t = [t; t5];

%% FINAL OUTPUTS

TotalT = sum(t);
TotalDV = sum(abs(dV));

fprintf('\nStandard strategy');
fprintf('\nTotal DV spent is: %f km/s', TotalDV);
fprintf('\nTime in flight: %f s\n', round(TotalT));

%% PLOT

OrbitsToPlot = [Orbit1; Orbit2; Orbit3; Orbit4; FinalOrbit];

orbitDraw(OrbitsToPlot)