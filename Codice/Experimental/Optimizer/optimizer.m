clear all; close all; clc;

%% DATA AND UTILS IMPORT 

addpath(genpath("..\Data\"))
load("PForbs.mat");
load("utils.mat",'mu');

%% DATA EXTRACTION

% Starting orbit
a1 = orb1(1);
e1 = orb1(2);
i1 = orb1(3);
O1 = orb1(4);
o1 = orb1(5);
r1p = a1*(1-e1);
r1a = a1*(1+e1);
p1 = a1*(1-e1^2);

% Final orbit
a2 = orb2(1);
e2 = orb2(2);
i2 = orb2(3);
O2 = orb2(4);
o2 = orb2(5);
r2p = a2*(1-e2);
r2a = a2*(1+e2);
p2 = a2*(1-e2^2);


%% FROM STARTING TO TRANSFER 1 ORBIT

% Optimal transfer orbit 1
syms rta
rtp = r1p;
et = (rta-rtp)/(rta+rtp);
at = (rta+rtp)/2;
pt = at * (1-et^2);

% Delta V
v1p = sqrt(mu/p1) * (1+e1);
vtp = sqrt(mu/pt) * (1+et);
dV1 = abs(vtp-v1p);

%% PERIAPSIS ALIGNMENT

% Transfer orbit 2
rt2p = r2p;
rt2a = rta;
at2 = (rt2a+rt2p)/2;
et2 = (rt2a-rt2p)/(rt2a+rt2p);
pt2 = at2 * (1-et2^2);

% Delta V
vta = sqrt(mu/pt) * (1+et*cos(pi));
vt2a = sqrt(mu/pt2) * (1+et2*cos(pi));
dV12 = abs(vta-vt2a);

%% PLANE CHANGE

% New orbit angles
dO = O2 - O1;
di = i2 - i1;

% Cases 
if dO > 0 && di > 0
    id = 1;
elseif dO > 0 && di < 0
    id = 2;
elseif dO < 0 && di > 0
    id = 3;
elseif dO < 0 && di < 0
    id = 4;
else 
    error("Error, do and/or di are not valid, please check inputs!")
end

% Spheric triangles
alpha = acos(cos(i1)*cos(i2) + sin(i1)*sin(i2)*cos(dO));

uiCos= (-1)^(id+1) * (cos(alpha)*cos(i1) - cos(i2))/(sin(alpha)*sin(i1));
uiSin = (-1)^(id+1) * sin(i2)*sin(dO)/sin(alpha);
ui = atan(uiSin/uiCos);

ufCos = (-1)^(id+1) * (cos(i1) - cos(alpha)*cos(i2))/(sin(alpha)*sin(i2));
ufSin = (-1)^(id+1) * sin(i1)*sin(dO)/sin(alpha);
uf = atan(ufSin/ufCos);


%% PLANE CHANGE
% Theta of intersection
thetai = wrapTo2Pi(ui - o1);

% Periapsis argument 
thetaf = thetai;
of = wrapTo2Pi(uf - thetaf);
% Delta V
Vt2t = sqrt(mu/pt2) * (1 + et2*cos(thetai));
dV2 = abs(2 * Vt2t * sin(alpha/2));

%% ARGUMENT OF PERIAPSIS ALIGNMENT

% Delta argument of periapsis 
do = o2-of;

% Delta V
vt2r = sqrt(mu/pt2)*et2*sin(do/2);
dV3 = abs(2*vt2r);

%% FROM TRANSFER 2 TO TARGET ORBIT

% Delta V
vt2p = sqrt(mu/pt2) * (1+et2);
v2p = sqrt(mu/p2) * (1+e2);
dV4 = abs(vt2p-v2p);

%% TOTAL IMPULSE

dv = dV1+dV2+dV3+dV12+dV4;

% From syms to function handler
dvFunc = matlabFunction(dv);
dV2Func = matlabFunction(dV2);

% Radii to consider
rs = r1a:1:100000;
dvs = arrayfun(dvFunc,rs);
[min,index] = min(dvs);
min
rs(index)

rta = rs(index);
dV2Func(rta)
orbits = [a1,e1,i1,O1,o1,orb1(6),2*pi;...
((rta+rtp)/2),(abs(rta-rtp)/(rta+rtp)),i1,O1,o1,0,pi;...
((rta+rt2p)/2),(abs(rta-rt2p)/(rta+rt2p)),i1,O1,o1,pi,thetai;...
((rta+rt2p)/2),(abs(rta-rt2p)/(rta+rt2p)),i2,O2,of,thetaf,(o2-of)/2;...
((rta+rt2p)/2),(abs(rta-rt2p)/(rta+rt2p)),i2,O2,o2,(o2-of)/2,2*pi;...
 a2,e2,i2,O2,o2,0,orb2(6)];

orbitDraw(orbits,0.01);
