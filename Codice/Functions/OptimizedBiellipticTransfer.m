clear 
close all 
clc

%% OPTIMIZED BIELLIPTICAL TRANSFER ORBIT 

% This matlab script calcultes an optimezed bielliptical transfer by giving
% it the initial and final orbits.
% It first reaches the pericenter and then shifts to a first transfer
% orbit that uses as a apocenter distance a variable cycling form the first
% obrit's radius to 100'000 Km. At the apocenter it moves to a second orbit
% which is tangent to the pericenter of the target point. After we
% perform a plane change and an argument of periapsis allignment on the
% second orbit and at the pericenter we move to the final orbit, where we
% wait to reach the final point. 
% The script is optimezed to calculate a transfer radius which corrisponds to
% the minimum velocity difference required to complete the transfer.

%% DATA AND UTILS IMPORT 

addpath(genpath("../Data/"))
load("GivenData.mat");
load("utils.mat",'mu');

%% DATA EXTRACTION

[a1,e1,i1,O1,o1,thi] = carToOrbital(r1,V1);

% Starting orbit
Orbit1(1) = a1;
Orbit1(2) = e1;
Orbit1(3) = i1;
Orbit1(4) = O1;
Orbit1(5) = o1;
Orbit1(6) = thi;
r1p = a1*(1-e1);
r1a = a1*(1+e1);
p1 = a1*(1-e1^2);

% Final orbit
a2 = FinalOrbit(1);
e2 = FinalOrbit(2);
i2 = FinalOrbit(3);
O2 = FinalOrbit(4);
o2 = FinalOrbit(5);
r2p = a2*(1-e2);
r2a = a2*(1+e2);
p2 = a2*(1-e2^2);
 

%% FROM STARTING POINT TO 1st TRANSFER ORBIT

Orbit1(7) = 0;

% Optimal first transfer orbit
syms rta
rtp = r1p;
et = (rta-rtp)/(rta+rtp);
at = (rta+rtp)/2;
pt = at * (1-et^2);

% First velocity difference
v1p = sqrt(mu/p1) * (1+e1);
vtp = sqrt(mu/pt) * (1+et);
dV1 = abs(vtp-v1p);

%% ALIGNMENT WITH TARGET ORBIT

% Second transfer orbit 
rt2p = r2p;
rt2a = rta;
at2 = (rt2a+rt2p)/2;
et2 = (rt2a-rt2p)/(rt2a+rt2p);
pt2 = at2 * (1-et2^2);

% Second velocity difference
vta = sqrt(mu/pt) * (1+et*cos(pi));
vt2a = sqrt(mu/pt2) * (1+et2*cos(pi));
dV2 = abs(vta-vt2a);

%% PLANE CHANGE

% Orbital Angles Difference
dO = O2 - O1;
di = i2 - i1;

% Study of possible outcomes
if dO > 0 && di > 0
    id = 1;
elseif dO > 0 && di < 0
    id = 2;
elseif dO < 0 && di > 0
    id = 3;
elseif dO < 0 && di < 0
    id = 4;
else 
    error("Error: either di or dO is null therefore the spheric angles equations can't be used. Please change inputs.")
end

% Spheric Angles Equations
alpha = acos(cos(i1)*cos(i2) + sin(i1)*sin(i2)*cos(dO));

u1Cos= (-1)^(id+1) * (cos(alpha)*cos(i1) - cos(i2))/(sin(alpha)*sin(i1));
u1Sin = (-1)^(id+1) * sin(i2)*sin(dO)/sin(alpha);
u1 = atan(u1Sin/u1Cos);

u2Cos = (-1)^(id+1) * (cos(i1) - cos(alpha)*cos(i2))/(sin(alpha)*sin(i2));
u2Sin = (-1)^(id+1) * sin(i1)*sin(dO)/sin(alpha);
u2 = atan(u2Sin/u2Cos);

% Intersection theta
theta3f = wrapTo2Pi(u1 - o1);

% Finding most convenient intersection (closest to periapsis)
if abs(theta3f-pi) > pi/2
    theta3f = theta3f + pi;
end

% Periapsis Argument of the new orbit 
theta4i = theta3f;
of = wrapTo2Pi(u2 - theta4i);

% Third velocity difference
Vt2t = sqrt(mu/pt2) * (1 + et2*cos(theta3f));
dV3 = abs(2 * Vt2t * sin(alpha/2));

%% ARGUMENT OF PERIAPSIS ALIGNMENT

% Delta argument of periapsis 
do = wrapTo2Pi(o2-of);

% Maneuver Point
if do/2 < pi
    if theta4i < do/2 || theta4i > do/2 + pi
        theta4f = do/2;
    else
        theta4f = do/2 + pi;
    end
else
    if theta4i < do/2 && theta4i > do/2 + pi
        theta4f = do/2;
    else
        theta4f = do/2 + pi;
    end
end

% New orbit's initial theta
if do/2 < pi
    if theta4i < do/2 || theta4i > do/2 + pi
        theta5i = 2*pi - do/2;
    else
        theta5i = pi - do/2;
    end
else
    if theta4i < do/2 && theta4i > do/2 + pi
        theta5i = 2*pi - do/2;
    else
        theta5i = pi - do/2;
    end
end

% Fourth velocity difference
vt2r = sqrt(mu/pt2)*et2*sin(do/2);
dV4 = abs(2*vt2r);

%% FROM LAST TRANSFER ORBIT TO TARGET ORBIT

% Fifth velocity difference
vt2p = sqrt(mu/pt2) * (1+et2);
v2p = sqrt(mu/p2) * (1+e2);
dV5 = abs(vt2p-v2p);

%% OUTPUTS

% Total Velocity Difference
dv = [dV1,dV2,dV3,dV4,dV5];
totalDV = sum(dv);

%% CYCLE FUNCTIONS

% From syms to function handler
dvFunc = matlabFunction(totalDV);

% Radii to consider
rs = r1a:1:300000;
dvs = arrayfun(dvFunc,rs);
[min,index] = min(dvs);
min;
rs(index)

rta = rs(index);
rt2a = rta;

%% FINAL RADIUS AND VELOCITY VECTORS

[r2,V2] = orbitalToCar(a2,e2,i2,O2,o2,FinalOrbit(7));

%% ORBIT VECTORS CREATION

% Second orbit vector creation
Orbit2 = zeros(1,7);
Orbit2(1) = ((rta+rtp)/2);
Orbit2(2) = (abs(rta-rtp)/(rta+rtp));
Orbit2(3) = i1;
Orbit2(4) = O1;
Orbit2(5) = o1;
Orbit2(6) = 0;
Orbit2(7) = pi;

% Third orbit vector creation
Orbit3 = zeros(1,7);
Orbit3(1) = ((rta+rt2p)/2);
Orbit3(2) = (abs(rta-rt2p)/(rta+rt2p));
Orbit3(3) = i1;
Orbit3(4) = O1;
Orbit3(5) = o1;
Orbit3(6) = pi;
Orbit3(7) = theta3f;

% Fourth orbit vector creation
Orbit4 = zeros(1,7);
Orbit4(1) = ((rta+rt2p)/2);
Orbit4(2) = (abs(rta-rt2p)/(rta+rt2p));
Orbit4(3) = i2;
Orbit4(4) = O2;
Orbit4(5) = of;
Orbit4(6) = theta4i;
Orbit4(7) = theta4f;


% Fifth orbit vector creation
Orbit5 = zeros(1,7);
Orbit5(1) = ((rta+rt2p)/2);
Orbit5(2) = (abs(rta-rt2p)/(rta+rt2p));
Orbit5(3) = i2;
Orbit5(4) = O2;
Orbit5(5) = o2;
Orbit5(6) = theta5i;
Orbit5(7) = 0;

%% TIME TAKEN TO COMPLETE MANEUVER

% Time taken to reach first orbit's pericenter
t1 = timeOfFlight(Orbit1,thi,0);
t = t1;

% Time taken to reach second orbit's apocenter
t2 = pi*sqrt(((rta+rtp)/2)^3/mu);
t = [t; t2];

% Time taken to reach third orbit's manuever point
t3 = timeOfFlight(Orbit3,pi,theta3f);
t = [t; t3];

% Time taken to reach fourth orbit's manuever point
t4 = timeOfFlight(Orbit4,Orbit4(1,6),Orbit4(1,7));
t = [t; t4];

% Time taken to reach fifth orbit's pericenter
t5 = timeOfFlight(Orbit5,Orbit5(1,6),Orbit5(1,7));
t = [t; t5];

% Time taken to reach final destination
FinalOrbit(1,6) = 0;
t6 = timeOfFlight(FinalOrbit,FinalOrbit(1,6),FinalOrbit(1,7));
t = [t; t6];

% Total Time Required
totalT = sum(t)

%% VELOCITIES DIFFERENCES

% First velocity difference
v1p = sqrt(mu/p1) * (1+e1);
vtp = sqrt(mu/(2*rta*rtp/(rta+rtp))) * (1+abs(rta-rtp)/(rta+rtp));
dV1 = abs(vtp-v1p);

% Second velocity difference
vta = sqrt(mu/(2*rta*rtp/(rta+rtp))) * (1+abs(rta-rtp)/(rta+rtp)*cos(pi));
vt2a = sqrt(mu/(2*rt2a*rt2p/(rt2a+rt2p))) * (1+abs(rta-rt2p)/(rta+rt2p)*cos(pi));
dV2 = abs(vta-vt2a);

% Third velocity difference
Vt2t = sqrt(mu/(2*rt2a*rt2p/(rt2a+rt2p))) * (1 + abs(rta-rt2p)/(rta+rt2p)*cos(theta3f));
dV3 = abs(2 * Vt2t * sin(alpha/2));

% Fourth velocity difference
vt2r = sqrt(mu/(2*rt2a*rt2p/(rt2a+rt2p)))*abs(rta-rt2p)/(rta+rt2p)*sin(do/2);
dV4 = abs(2*vt2r);

% Fifth velocity difference
vt2p = sqrt(mu/(2*rt2a*rt2p/(rt2a+rt2p))) * (1+abs(rta-rt2p)/(rta+rt2p));
v2p = sqrt(mu/p2) * (1+e2);
dV5 = abs(vt2p-v2p);

% Total Velocity Difference
dv = [dV1,dV2,dV3,dV4,dV5];
totalDV = sum(dv)

%% PLOT

orbits = [Orbit1; Orbit2; Orbit3; Orbit4; Orbit5; FinalOrbit];


% drawOrbitapp(orbits,0.01,1,[54,13],'..\..\Images\Optimized\optimized')
% drawOrbitapp(orbits,0.01,1,[105,13],'..\..\Images\Optimized\optimized2')
orbitDraw(orbits);

