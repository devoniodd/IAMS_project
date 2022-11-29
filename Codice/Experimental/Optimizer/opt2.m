clear all; 
close all; clc;

%% DATA AND UTILS IMPORT 

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

%% CIRCULARIZATION

v1a = sqrt(mu/p1) * (1-e1);
v1c = sqrt(mu/r1a);

dV1 = abs(v1c-v1a);

%% TO HIGHER ORBIT

% Optimal transfer orbit 1
syms rta
rtp = r1a;
et = abs(rta-rtp)/(rta+rtp);
at = (rta+rtp)/2;
pt = at * (1-et^2);

% Delta V
vtp = sqrt(mu/pt) * (1+et);
dV2 = abs(vtp-v1c);

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

% Delta omega
do = (ui-pi)-o1;
oc = o1 + do;

% Theta of intersection
thetai = wrapTo2Pi(ui - oc);

% Periapsis argument 
thetaf = thetai;
of = wrapTo2Pi(uf - thetaf);

% Delta V
Vtt = sqrt(mu/pt) * (1 + et*cos(thetai));
dV3 = abs(2 * Vtt * sin(alpha/2));

%% LOWER 

% New orbit
rt2p = r2a;
rt2a = rta;
at2 = (rt2a+rt2p)/2;
et2 = abs(rt2a-rt2p)/(rt2a+rt2p);
pt2 = at2 * (1-et2^2);

% Adapt apoapsis
vta = sqrt(mu/pt) * (1-et);
vt2a = sqrt(mu/pt2) * (1-et2);
dV4 = abs(vt2a-vta);

% Cicrcular
v2c = sqrt(mu/r2a);
vt2p =  sqrt(mu/pt2) * (1+et2);
dV5 = abs(v2c-vt2p);

%% ARGUMENT OF PERIAPSIS ALIGNMENT

% Delta argument of periapsis 
do = o2-of;

%% SHAPE ADAPT

% Adapt periapsis
v2a = sqrt(mu/p2) * (1-e2);
dV6 = abs(v2c-v2a);

%% OPT VELOCITY

dV = dV1+dV2+dV3+dV4+dV5+dV6;

% From syms to function handler
dvFunc = matlabFunction(dV);

% Radii to consider
rs = r1a:10:1000000;
dvs = arrayfun(dvFunc,rs);
[min,index] = min(dvs);
min
rs(index)

%% PLOTS

figure;
fplot(dvFunc,[r1p 100000],'lineWidth',2);
xline(r1a);
xline(r2a);
grid on;
legend("Total dV","Orbit 1 apogee radius","Orbit 2 apogee radius",Location="best");
ylim('auto');
hold off;

dV2Func = matlabFunction(dV2);
dV3Func = matlabFunction(dV3);
dV4Func = matlabFunction(dV4);
dV5Func = matlabFunction(dV5);

figure;
fplot(dV2Func,[r2a 600000],'lineWidth',2);
hold on;
fplot(dV3Func,[r2a 600000],'lineWidth',2);
fplot(dV4Func,[r2a 600000],'lineWidth',2);
fplot(dV5Func,[r2a 600000],'lineWidth',2);
legend("To first elliptic orbit dV","Plane change dV","To second elliptic orbit dV","Final orbit circularization dV",Location="best");

