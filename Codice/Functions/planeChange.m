function [dV,tof,targetOrbit,currentOrbit] = planeChange(currentOrbit,i2,O2)

%% UTILS
if ismac
    load("../Data/utils.mat",'mu');
else
    load("..\Data\utils.mat",'mu');
end

%% INPUT ORGANIZATION
a1 = currentOrbit(1);
e1 = currentOrbit(2);
i1 = currentOrbit(3);
O1 = currentOrbit(4);
o1 = currentOrbit(5);

%% PROBLEM VARIABLES
dO = O2 - O1;
di = i2 - i1;

%% CASES 

if dO > 0 && di > 0
    id = 1;
elseif dO > 0 && di < 0
    id = 2;
elseif dO < 0 && di > 0
    id = 3;
elseif dO < 0 && di < 0
    id = 4;
else 
    error("Error, i2 and/or O2 are not valid, please check inputs!")
end

%% SPHERIC TRIANGLES
alpha = acos(cos(i1)*cos(i2) + sin(i1)*sin(i2)*cos(dO));

u1Cos= (-1)^(id+1) * (cos(alpha)*cos(i1) - cos(i2))/(sin(alpha)*sin(i1));
u1Sin = (-1)^(id+1) * sin(i2)*sin(dO)/sin(alpha);
u1 = atan(u1Sin/u1Cos);

u2Cos = (-1)^(id+1) * (cos(i1) - cos(alpha)*cos(i2))/(sin(alpha)*sin(i2));
u2Sin = (-1)^(id+1) * sin(i1)*sin(dO)/sin(alpha);
u2 = atan(u2Sin/u2Cos);

theta1 = wrapTo2Pi(u1-o1);

theta2 = theta1;
o2 =wrapTo2Pi(u2 - theta2);

%% dV
V1t = sqrt(mu/(a1*(1-e1^2))) * (1 + e1*cos(theta1));
dV = 2 * V1t * sin(alpha/2);

%% OUTPUT

% Current Orbit
currentOrbit(1,7) = theta1;

% Target Orbit
targetOrbit = zeros(1,7);
targetOrbit(1,1) = a1;
targetOrbit(1,2) = e1;
targetOrbit(1,3) = i2;
targetOrbit(1,4) = O2;
targetOrbit(1,5) = o2;
targetOrbit(1,6) = theta2;
targetOrbit(1,7) = nan;

% Time to reach maneuver point
tof = timeOfFlight(currentOrbit,currentOrbit(1,6),currentOrbit(1,7));








