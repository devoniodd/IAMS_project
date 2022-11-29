clear all;
close all;
clc;

addpath(genpath("../../Data/"))
load("utils.mat",'mu');

r1 = [-8048.2861, -4171.3048, 2895.9296];
v1 = [1.7540, -5.9910, -1.9520];

orbit1 = zeros(1,7);
[orbit1(1), orbit1(2), orbit1(3), orbit1(4), orbit1(5), orbit1(6)] = carToOrbital(r1,v1);
e1 = orbit1(2);
a1 = orbit1(1);
p1 = a1 * (1-e1^2);
theta1 = orbit1(6);
i1 = orbit1(3);
O1 = orbit1(4);
v1t = sqrt(mu/p1) * (1+e1*cos(theta1));
v1r = sqrt(mu/p1) * (e1*sin(theta1));

orbit2 = [16720.0000 0.2502 1.1190 0.6245 3.1350 3.1000];
[r2, v2] = orbitalToCar(16720.0000,0.2502,1.1190,0.6245,3.1350,3.1000);
r2 = r2';

t = 4000;
a = 4000;
e = 0.5;
vLim = 7;
orbitT = zeros(1,7);
dVlambert = vLim;
while a * (1-e) < 6578
t = t+10;
vt1= lambert(r1,r2,t,1);
[a, e, i, O] = carToOrbital(r1,vt1);
end
tts = [];
vvs = [];
while dVlambert >= vLim
    t = t+100;
    [vt1, vt2] = lambert(r1,r2,t,1);
    [at, et, it, Ot, ~, thetat] = carToOrbital(r1,vt1);
    pt = at * (1-et^2);
    vtt = sqrt(mu/pt) * (1+et*cos(thetat));
    vtr = sqrt(mu/pt) * (et*sin(thetat));
    alpha = acos(cos(i1)*cos(it) + sin(i1)*sin(it)*cos(O1-Ot));
    dVlambert = sqrt((vtr-v1r)^2+v1t^2+vtt^2-2*vtt*v1t*cos(alpha));
    tts = [tts t];
    vvs = [vvs dVlambert];
end

[orbitT(1), orbitT(2), orbitT(3), orbitT(4), orbitT(5), orbitT(6)] = carToOrbital(r1,vt1);
[~, ~, ~, ~, ~, orbitT(7)] = carToOrbital(r2,vt2);
orbit1(7) = 2*pi;
orbit2(7) = 2*pi;
orbit2(6) = 0;
orbits = [orbit1; orbitT; orbit2];

%% Orbits

Ot = orbitT(4);
ot = orbitT(5);
it = orbitT(3);
et = orbitT(2);
at = orbitT(1);
pt = at * (1-et^2);
thetat = orbitT(6);

O2 = orbit2(4);
ot = orbit2(5);
i2 = orbit2(3);
e2 = orbit2(2);

%% dV Lambert 


dO = O2 - Ot;
di = i2 - it;

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
alpha = acos(cos(it)*cos(i2) + sin(it)*sin(i2)*cos(dO));

uiCos= (-1)^(id+1) * (cos(alpha)*cos(it) - cos(i2))/(sin(alpha)*sin(it));
uiSin = (-1)^(id+1) * sin(i2)*sin(dO)/sin(alpha);
ui = atan(uiSin/uiCos);

ufCos = (-1)^(id+1) * (cos(it) - cos(alpha)*cos(i2))/(sin(alpha)*sin(i2));
ufSin = (-1)^(id+1) * sin(it)*sin(dO)/sin(alpha);
uf = atan(ufSin/ufCos);

% Theta of intersection
thetai = wrapTo2Pi(ui - ot);

% Find most convenient intersection
if abs(thetai-pi) > pi/2
    thetai = thetai + pi;
end

% Periapsis argument 
thetaf = thetai;
of = wrapTo2Pi(uf - thetaf);

% Delta V
Vt2t = sqrt(mu/pt) * (1 + et*cos(thetai));
dV2 = abs(2 * Vt2t * sin(alpha/2));



%% ALIGN APSE LINE

orbitT2 = orbitT;
orbitT2(3) = orbit2(3);
orbitT2(4) = orbit2(4);
orbitT2(5) = of;

do = orbitT2(5) - orbit2(5);
[dVSecant, tof] = secantManeuver(orbitT2,orbit2,do,1);

dVlambert + dVSecant + dV2

t+tof


