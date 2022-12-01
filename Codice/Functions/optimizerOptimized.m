clear  
close all 
clc
%% OPTIMIZED STERNFELD MANOUVER

% This matlab script calcultes an optimezed bi-elliptical transfer by giving
% it the initial and final orbits.
% It first circularizes the starting orbit in order to chnage its argument
% of periapsis freely. It then moves to a second orbit in which the apogee
% radius is cycled between the apoapsis of the first orbit and 10'000'000 km.
% When the object reaches this said point it performs a plane change and
% switches to a second transfer orbit which has as the pericenter radius
% the apocenter radius of the final orbit. Once we arrive at the pericenter
% of this third orbit we perform a second circularizartion until we reach
% the tangent point with the target orbit. Once we reach this point we move
% onto the final orbit and wait until we reach our destination.
%

%% PLOT 

plot = 0;

%% DATA AND UTILS IMPORT 

addpath(genpath("../Data"));
load("PForbs.mat");
load("utils.mat",'mu');

%% DATA EXTRACTION

% Starting orbit
a1 = orb1(1);
e1 = orb1(2);
i1 = orb1(3);
O1 = orb1(4);
o1 = orb1(5);
theta1 = orb1(6);
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

%% FIRST CIRCULARIZATION

v1a = sqrt(mu/p1) * (1-e1);
v1c = sqrt(mu/r1a);

dV1 = abs(v1c-v1a);


%% FIRST TRANSFER ORBIT

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
ui = wrapTo2Pi(atan(uiSin/uiCos));

ufCos = (-1)^(id+1) * (cos(i1) - cos(alpha)*cos(i2))/(sin(alpha)*sin(i2));
ufSin = (-1)^(id+1) * sin(i1)*sin(dO)/sin(alpha);
uf = wrapTo2Pi(atan(ufSin/ufCos));

% Delta omega
do1 = wrapTo2Pi((ui-pi)-o1);
oc = wrapTo2Pi(o1 + do1);

% Theta of intersection
thetai = wrapTo2Pi(ui - oc);

% Periapsis argument 
thetaf = thetai;
of = wrapTo2Pi(uf - thetaf);

%% SECOND TRANSFER ORBIT 

% New orbit
rt2p = r2a;
rt2a = rta;
at2 = (rt2a+rt2p)/2;
et2 = abs(rt2a-rt2p)/(rt2a+rt2p);
pt2 = at2 * (1-et2^2);

% Delta V
Vtt = sqrt(mu/pt) * (1 + et*cos(thetai));
Vtr = sqrt(mu/pt) * et * sin(thetai);

Vt2t = sqrt(mu/pt2) * (1 + et2*cos(thetai));
Vt2r = sqrt(mu/pt2) * et2 * sin(thetai);

dV3 = sqrt((Vt2r-Vtr)^2 + Vt2t^2 + Vtt^2 - 2*Vtt*Vt2t*cos(alpha) );

%% SECOND CIRCULARIZATION

v2c = sqrt(mu/r2a);
vt2p =  sqrt(mu/pt2) * (1+et2);
dV4 = abs(v2c-vt2p);

%% ARGUMENT OF PERIAPSIS ALIGNMENT

% Delta argument of periapsis 
do = o2-of;

%% TARGET ORBIT ADAPTATION

% Adapt periapsis
v2a = sqrt(mu/p2) * (1-e2);
dV5 = abs(v2c-v2a);

%% OPTIMAL VELOCITY

dV = dV1+dV2+dV3+dV4+dV5;

% From syms to function handler
dvFunc = matlabFunction(dV);

% Radii to consider
[rOpt, dVmin] = fminbnd(dvFunc,r1a,10000000);

%% ORBITS' VECTORS CREATION

orb1(7) = pi;
orb1c1 = [r1a,0,i1,O1,o1,pi,do1];
orbT = [(rOpt+r1a)/2, abs(rOpt-r1a)/(rOpt+r1a), i1, O1, oc, 2*pi, thetai];
orbT2 = [(rOpt+r2a)/2, abs(rOpt-r2a)/(rOpt+r2a), i2, O2, of, thetai, 2*pi];
orb2c2 = [r2a, 0, i2, O2, o2,wrapTo2Pi(-do), pi];
orb2(6) = pi;

orbits = [orb1;orb1c1;orbT;orbT2;orb2c2;orb2];

%% DELTA Vs

dV2Func = matlabFunction(dV2);
dV3Func = matlabFunction(dV3);
dV4Func = matlabFunction(dV4);

dV = [dV1,dV2Func(rOpt),dV3Func(rOpt),dV4Func(rOpt),dV5];

%% TIME TAKEN TO COMPLETE MANEUVER

% Time taken to reach first orbit's apocenter
t1 = timeOfFlight(orb1,orb1(6),pi);
t = t1;

% Time taken to reach second orbit's apocenter
t2 = timeOfFlight(orb1c1,pi,orb1c1(7));
t = [t; t2];

% Time taken to reach third orbit's manuever point
t3 = timeOfFlight(orbT,2*pi,thetai);
t = [t; t3];

% Time taken to reach fourth orbit's manuever point
t4 = timeOfFlight(orbT2,thetai,2*pi);
t = [t; t4];

% Time taken to reach fifth orbit's pericenter
t5 = timeOfFlight(orb2c2,orb2c2(6),pi);
t = [t; t5];

% Time taken to reach final destination
t6 = timeOfFlight(orb2,pi,orb2(7));
t = [t; t6];

TotalT = sum(t);

%% PLOTS

if plot

    f1 = figure();
    fplot(dvFunc,[r1p 100000],'lineWidth',2);
    xl = xline(rOpt,'--',"$r_{opt}$",'Interpreter','latex');
    xl.FontSize = 18;
    xl.LineWidth = 1;
    xl.LabelVerticalAlignment = "top";
    xl.LabelOrientation = 'horizontal';
    grid on;
    legend("Total dV",'Location',"best",'Interpreter','latex','FontSize',15);
    ylim([4.5 7]);
    xticks([1e4,rOpt,3e4:1e4:10e4]);
    hold off;
    xlabel("\textbf{first transfer orbit apogee radius}",'Interpreter','latex','FontSize',15)
    ylabel("\textbf{total dV}",'Interpreter','latex','FontSize',15)

    f2 = figure();
    fplot(dV2Func,[r1a 100000],'lineWidth',2);
    hold on;
    fplot(dV3Func,[r1a 100000],'lineWidth',2);
    fplot(dV4Func,[r1a 100000],'lineWidth',2);
    xl = xline(rOpt,'--',"$r_{opt}$",'Interpreter','latex');
    xl.FontSize = 18;
    xl.LineWidth = 1;
    xl.LabelVerticalAlignment = "top";
    xl.LabelOrientation = 'horizontal';
    xticks([1e4,rOpt,3e4:1e4:10e4]);
    grid on;
    xlabel("\textbf{first transfer orbit apogee radius}",'Interpreter','latex','FontSize',15)
    ylabel("\textbf{dV}",'Interpreter','latex','FontSize',15)
    legend("to $1^{st}$ transfer orbit","plane and shape change","$2^{nd}$ tranfer orbit circularization",'Location',"northeast",'Interpreter','latex','FontSize',12.5);

    % Save
    saveas(f1,"../../Images/OptimizedOpt/totaldV",'png');
    saveas(f2,"../../Images/OptimizedOpt/partialdVs",'png');

    orbitDraw(orbits,[71 11],'../../Images/OptimizedOpt/optimizedOpt');

end

