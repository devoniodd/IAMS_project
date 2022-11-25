%% FAST STRATEGY 1
clear
close all
clc

% Manouver is exectuted with:
% Change of plane in the first aveilable manouver node
% Change of shape throug secant orbit
% Adaptation to final orbit

O_start =    [9723.68854055203  0.0865268303717344  0.461525983980154  1.17569878936122  0.982554186272932  1.40668276722613  6.28318530717959];
O_end =      [16720  0.250200000000000  1.11900000000000  0.624500000000000  3.13500000000000  0  pi];



%% FAST MANOUVER 2
clear
close all
clc

O_start =    [9723.68854055203  0.0865268303717344  0.461525983980154  1.17569878936122  0.982554186272932  1.40668276722613  6.28318530717959];
O_end =      [16720  0.250200000000000  1.11900000000000  0.624500000000000  3.13500000000000  0  pi];

% In order to save time the first orbit reduces its shape 
[dV,tof,currentOrbit,targetOrbit] = planeChange(O_start,O_end(3),O_end(4));
O_start(7)= O_start(6);
targetOrbit(7) = 2*pi;

rm0 = orbitalToCar(O_start(1),O_start(2),O_start(3),O_start(4),O_start(5),O_start(6));

ri = orbitalToCar(O_start(1),O_start(2),O_start(3),O_start(4),O_start(5),O_start(6));
rm1 = orbitalToCar(O_start(1),O_start(2),O_start(3),O_start(4),O_start(5),currentOrbit(7));
rm1 = rm1/norm(rm1) * 7000;

Orbit1 = orbitfrom2points(ri,pi,rm1);

% Executing plane and shape change
rm2 = orbitalToCar(O_end(1),O_end(2),O_end(3),O_end(4),O_end(5),2*pi);
rpend = orbitalToCar(O_end(1),O_end(2),O_end(3),O_end(4),O_end(5),0);
thetaM3 = 2*pi - acos(dot(rm0,rpend)/(norm(rm0)*norm(rpend)));
rm3 = orbitalToCar(O_end(1),O_end(2),O_end(3),O_end(4),O_end(5),thetaM3);
Orbit2 = orbitfrom2points(rm1,0,rm3);
O_end(6) = thetaM3;


Orbits = [O_start;Orbit1;Orbit2;O_end];

% DV CALCULATION
DV1 = simpleDVcalculator(O_start,Orbit1);
DV2 = simpleDVcalculator(Orbit1,Orbit2);
DV3 = simpleDVcalculator(Orbit2,O_end);

DV_tot = DV1+DV2+DV3;

% TIME CALCULATION
Dt1 = timeOfFlight(Orbit1);
Dt2 = timeOfFlight(Orbit2);
Dt3 = timeOfFlight(O_end);

Time = Dt1 + Dt2 + Dt3;

fprintf('\nTime optimized manouver 2');
fprintf('\nTotal DV spent is: %f km/s', DV_tot);
fprintf('\nTime in flight: %f s\n', round(Time));

% PLOT
drawOrbitapp(Orbits,0.1);
[lat,lon] = orbitpropagator(Orbits,15,10000);