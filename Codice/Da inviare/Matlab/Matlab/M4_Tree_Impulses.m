delete(gcp('nocreate')); % Close every multi core environment
clear
close all
clc

%% 3 IMPULSE MANOUVER
% DESCRIPTION TBA

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

% CHANGE SEMI MAJOR APSIS Orbit 2
[dV,tof,hypoteticalorbit,targetorbit] = planeChange(O_end,O_start(3),O_start(4),2); % Finding theta on final orbit
thetaf = hypoteticalorbit(7);

[dV,tof,Orbit1,targetorbit] = planeChange(Orbit1,O_end(3),O_end(4)); % Finding manouver node for first manouver
rAt = O_end(1)*(1-O_end(2)^2)/(1+O_end(2)*cos(thetaf));              % Apogee of transfer orbit
a2 = (Orbit1(1) + rAt)/2;
e2 = (- Orbit1(1) + rAt)/(2*a2);
omega2 = pi + Orbit1(5) + Orbit1(7) - Orbit1(6);
Orbit2 = [a2,e2,Orbit1(3),Orbit1(4),omega2,0,pi];
O_end(6) = thetaf;

% DV CALCULATION
DV = zeros(1,3);
DV(1) = simpleDVcalculator(O_start,Orbit1);
DV(2) = simpleDVcalculator(Orbit1,Orbit2);
DV(3) = simpleDVcalculator(Orbit2,O_end);

DV_tot = sum(DV);

% TIME CALCULATION
Dt = zeros(1,4);
Dt(1) = timeOfFlight(O_start);
Dt(2) = timeOfFlight(Orbit1);
Dt(3) = timeOfFlight(Orbit2);
Dt(4) = timeOfFlight(O_end);

Time = sum(Dt);

fprintf('\nTime optimized manouver 2');
fprintf('\nTotal DV spent is: %f km/s', DV_tot);
fprintf('\nTime in flight: %f s\n', round(Time));

orbits = [O_start;Orbit1;Orbit2;O_end];
orbitDraw(orbits);

%% DISPLAYING DV VECTORS AT 3rd MANOUVER NODE
% SETTING UNIT VECTORS FOR DIRECTIONS
[rm3_1,v3_1] = orbitalToCar(Orbit2(1),Orbit2(2),Orbit2(3),Orbit2(4),Orbit2(5),Orbit2(7));
[rm3_2,v3_2] = orbitalToCar(O_end(1),O_end(2),O_end(3),O_end(4),O_end(5),O_end(6));

rDir1 = - rm3_1/norm(rm3_1);
rDir2 = - rm3_2/norm(rm3_2);
hDir1 = cross(rDir1,v3_1);
hDir1 = hDir1/norm(hDir1);
hDir2 = cross(rDir2,v3_2);
hDir2 = hDir2/norm(hDir2);

thDir1 = cross(rDir1,hDir1);
thDir1 = - thDir1/norm(thDir1);
thDir2 = cross(rDir2,hDir2);
thDir2 = - thDir2/norm(thDir2);

i = [1,0,0];
j = [0,1,0];
k = [0,0,1];

% CALCULATING DVs
% Importing variables
load("..\Data\utils.mat",'mu');

% Calculating velocity
p1 = Orbit2(1)*(1 - Orbit2(2)^2);
p2 = O_end(1)*(1 - O_end(2)^2);

Vtheta1 = sqrt(mu/p1)*(1 - Orbit2(2));
Vr1 = 0;    % Manouver is executed at apogee
Vtheta2 = sqrt(mu/p2)*(1 + O_end(2)*cos(O_end(6)));
Vr2 = sqrt(mu/p2)*O_end(2)*sin(O_end(6));

% Calculating planechange DV
gamma = acos(cos(Orbit2(3)) * cos(O_end(3)) + sin(Orbit2(3))*sin(O_end(3))*cos(O_end(4)-Orbit2(4)));
DVchangeplane = - thDir1 * Vtheta1*(1 - cos(gamma)) + hDir1 * (Vtheta1*sin(gamma));

% Calculating change shape and periapsis DV
DVtheta = thDir2 * (Vtheta2 - Vtheta1);
DVrad = - rDir2 * (Vr2);

DVreal = (v3_2-v3_1);

% SETTING PLOT ENVIRONMENT
figure Name 'Speed Vector Plot'
hold on;
axis equal;
grid on;

xlim([-2,2]);
ylim([-2,2]);
zlim([-5,2]);

DVtotvett = DVchangeplane + DVtheta + DVrad;

set(gca,'color','k','xcolor','w','ycolor','w','zcolor','w');
set(gcf,'color','k');

%% ADDING ORBIT
color = jet(4);

itmax = 100;
thetaiOrb = linspace(Orbit2(7)-0.1,Orbit2(7),itmax);
thetaiComp = linspace(Orbit2(7),Orbit2(7)+0.1,itmax);
orb2 = zeros(3,itmax);
orb2comp = zeros(3,itmax);

thetafOrb = linspace(O_end(6),O_end(6)+0.1,itmax);
thetafComp = linspace(O_end(6)-0.1,O_end(6),itmax);
orb3 = zeros(3,itmax);
orb3comp = zeros(3,itmax);

for n = 1: itmax
    orb2(:,n) = orbitalToCar(Orbit2(1),Orbit2(2),Orbit2(3),Orbit2(4),Orbit2(5),thetaiOrb(n));
    orb2(:,n) = orb2(:,n) - rm3_2;
    orb2comp(:,n) = orbitalToCar(Orbit2(1),Orbit2(2),Orbit2(3),Orbit2(4),Orbit2(5),thetaiComp(n));
    orb2comp(:,n) = orb2comp(:,n) - rm3_2;
    orb3(:,n) = orbitalToCar(O_end(1),O_end(2),O_end(3),O_end(4),O_end(5),thetafOrb(n));
    orb3(:,n) = orb3(:,n) - rm3_1;
    orb3comp(:,n) = orbitalToCar(O_end(1),O_end(2),O_end(3),O_end(4),O_end(5),thetafComp(n));
    orb3comp(:,n) = orb3comp(:,n) - rm3_1;
end

% line([0,DVreal(1)],[0,DVreal(2)],[0,DVreal(3)],'Color','y','linewidth', 2, 'linestyle','-.')
h = plot3(nan,nan,nan,'.',MarkerSize=22,Color='w');
set(h,'XData',0,'YData',0,'ZData',0);

Orbit2 = plot3(orb2(1,:),orb2(2,:),orb2(3,:),'LineWidth' , 1.5, 'Color', color(3,:));
plot3(orb2comp(1,:),orb2comp(2,:),orb2comp(3,:),'LineStyle','--', 'Color', color(3,:));
Orbitend = plot3(orb3(1,:),orb3(2,:),orb3(3,:),'LineWidth' , 1.5, 'Color', color(4,:));
plot3(orb3comp(1,:),orb3comp(2,:),orb3comp(3,:),'LineStyle','--', 'Color', color(4,:));

%% ADDING VECTORS
DVCP = line([0,DVchangeplane(1)],[0,DVchangeplane(2)],[0,DVchangeplane(3)],'Color','cyan','linewidth', 2);
DVCS = line([0,DVtheta(1)],[0,DVtheta(2)],[0,DVtheta(3)],'Color','b','linewidth', 2);
DVCS2 = line([0,DVrad(1)],[0,DVrad(2)],[0,DVrad(3)],'Color','g','linewidth', 2);
DVTOT =line([0,DVtotvett(1)],[0,DVtotvett(2)],[0,DVtotvett(3)],'Color','w','linewidth', 2);
labels = legend([DVCP,DVCS,DVCS2,DVTOT,Orbit2,Orbitend],'Change Plane $\Delta$V','Shape/periapsis change tangent $\Delta$V','Shape/periapsis change radial $\Delta$V','Total $\Delta$V','Second Orbit','Final Orbit');
labels.Interpreter = 'latex';
labels.TextColor = 'w';
labels.FontSize = 15;
