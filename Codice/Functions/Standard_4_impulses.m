delete(gcp('nocreate')); % Close every multi core environment
clear
close all
clc

%% 3 IMPULSE MANOUVER

% This matlab script calcultes an optimized 3 impulse manouver with
% relative DV spent and time of flight.
% The first manouver consist of a circularization at apoapsis, the
% second is a DV optimized direct transfer from the intersection between
% the two orbital planes to a generic node located on the final orbit.
% Please note that the second manouver will change more than 1 orbital
% parameter at a time, therefore for ease of coding (not necesserely for 
% optimization), transfer orbits will be calculated with the function 
% Orbitfrom2point.m, and DV will be obtained vectorially with the function 
% simpleDVcalculator.m 

%% DATA AND UTILS IMPORT

addpath(genpath("../Data/"))
load("GivenData.mat");
load("utils.mat",'mu');

[aS,eS,iS,OmegaS,omegaS,thetaS] = carToOrbital(r1,V1);
O_start(1) = aS;
O_start(2) = eS;
O_start(3) = iS;
O_start(4) = OmegaS;
O_start(5) = omegaS;
O_start(6) = thetaS;
O_start(7) = nan;

O_end = FinalOrbit;

% CIRCULARIZE AT APOAPSIS Orbit 1
O_start(7) = pi;
a0 = O_start(1);
e0 = O_start(2);
rA0 = a0*(1+e0);

Orbit1 = O_start;
Orbit1(1) = rA0;
Orbit1(2) = 0;
Orbit1(6) = pi;

% CHANGE SHAPE Orbit 2
[dV,tof,Orbit1,targetorbit] = planeChange(Orbit1,O_end(3),O_end(4)); % Finding manouver node
rPend = O_end(1)*(1-O_end(2));
a2 = (Orbit1(1) + rPend)/2;
e2 = (- Orbit1(1) + rPend)/(2*a2);
omega2 = pi + Orbit1(5) + Orbit1(7) - Orbit1(6);
Orbit2 = [a2,e2,Orbit1(3),Orbit1(4),omega2,0,pi];

% CHANGE OF PLANE Orbit 3
[dV,tof,Orbit2,Orbit3] = planeChange(Orbit2,O_end(3),O_end(4),2); % Executing manouver;
Orbit3(7) = pi;

% CHANGE OF PLANE SEMI-MAJOR APSIS Orbit 4
[dV,Orbit3,Orbit4] = changeSemiMajorApsis(Orbit3,O_end(1));

% CHANGE OF PERIAPSIS
[dV,tof,Orbit4,Orbit5] = changePeriapsisArg(Orbit4,O_end(5));
Orbit5(7) = O_end(7);

% DV CALCULATION
DV = zeros(1,4);
DV(1) = simpleDVcalculator(O_start,Orbit1);
DV(2) = simpleDVcalculator(Orbit1,Orbit2);
DV(3) = simpleDVcalculator(Orbit2,Orbit4);
DV(4) = simpleDVcalculator(Orbit4,Orbit5);

DV_tot = sum(DV);

% TIME CALCULATION
Dt = zeros(1,5);
Dt(1) = timeOfFlight(O_start);
Dt(2) = timeOfFlight(Orbit1);
Dt(3) = timeOfFlight(Orbit2);
Dt(4) = timeOfFlight(Orbit4);
Dt(5) = timeOfFlight(Orbit5);

Time = sum(Dt);

fprintf('\nTime optimized manouver 2');
fprintf('\nTotal DV spent is: %f km/s', DV_tot);
fprintf('\nTime in flight: %f s\n', round(Time));

orbits = [O_start;Orbit1;Orbit2;Orbit4;Orbit5];
orbitDraw(orbits);

