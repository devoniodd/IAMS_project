%% FAST STRATEGY 1
clear
close all
clc

% Manouver is exectuted with:
% Change of plane in the first aveilable manouver node
% Change of shape throug secant orbit
% Adaptation to final orbit

O_start =    [9723.68854055203  0.0865268303717344  0.461525983980154  1.17569878936122  0.982554186272932  0 2*pi];
O_end =      [16720  0.250200000000000  1.11900000000000  0.624500000000000  3.13500000000000  0  2*pi];

omega_start = O_start(4);

% FINDING DIRECTION OF H VECTOR OF INITIAL ORBIT
thstart = O_start(3);
phistart = wrapTo2Pi(O_start(4)-pi/2);
hDir1 = [cos(phistart)*sin(thstart),sin(phistart)*sin(thstart),cos(thstart)];
hDir1 = hDir1/norm(hDir1);

% FINDING DIRECTION OF H VECTOR OF TARGET ORBIT
thend = O_end(3);
phiend = wrapTo2Pi(O_end(4)-pi/2);
hDir2 = [cos(phiend)*sin(thend),sin(phiend)*sin(thend),cos(thend)];
hDir2 = hDir2/norm(hDir2);

% DIRECTION OF NODAL APSIS
nodAps = cross(hDir2,hDir1);

%% FAST MANOUVER 2 - Reduce a before manouver
% MANOUVER LIST%
% [CHANGE OF SHAPE + PERIAXIS]
% [CHANGE OF PLANE + CHANGE OF SHAPE]
% [CHANGE OF SHAPE + ]
clear
close all
clc

O_start =    [9723.68854055203  0.0865268303717344  0.461525983980154  1.17569878936122  0.982554186272932  1.40668276722613  6.28318530717959];
O_end =      [16720  0.250200000000000  1.11900000000000  0.624500000000000  3.13500000000000  0  3.1];

% In order to save time the first orbit reduces its shape 
[dV,tof,currentOrbit,targetOrbit] = planeChange(O_start,O_end(3),O_end(4));
O_start(7)= O_start(6);
targetOrbit(7) = 2*pi;

rm0 = orbitalToCar(O_start(1),O_start(2),O_start(3),O_start(4),O_start(5),O_start(6));

ri = orbitalToCar(O_start(1),O_start(2),O_start(3),O_start(4),O_start(5),O_start(6));
rm1 = orbitalToCar(O_start(1),O_start(2),O_start(3),O_start(4),O_start(5),currentOrbit(7));
rm1 = rm1/norm(rm1) * 14000;

Orbit1 = orbitfrom2points(ri,rm1,"theta0",0);

% Executing plane and shape change
rpend = orbitalToCar(O_end(1),O_end(2),O_end(3),O_end(4),O_end(5),0);
thetaM3 = 2*pi - acos(dot(rm0,rpend)/(norm(rm0)*norm(rpend)));
rm3 = orbitalToCar(O_end(1),O_end(2),O_end(3),O_end(4),O_end(5),thetaM3);
Orbit2 = orbitfrom2points(rm1,rm3,"theta0",pi);
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
orbitDraw(Orbits);
%orbitpropagator(Orbits,15,1000,[1,1,0],"peri");

%% FAST MANOUVER 2.1.2 - Apoapsis at node a before manouver
% MANOUVER LIST%
% [CHANGE OF SHAPE + PERIAXIS]
% [CHANGE OF PLANE + CHANGE OF SHAPE]
% [CHANGE OF SHAPE + ]
clear
close all
clc

O_start =    [9723.68854055203  0.0865268303717344  0.461525983980154  1.17569878936122  0.982554186272932  1.40668276722613  6.28318530717959];
O_end =      [16720  0.250200000000000  1.11900000000000  0.624500000000000  3.13500000000000  0  3.1];

% In order to save time the first orbit reduces its shape 
[dV,tof,currentOrbit,targetOrbit] = planeChange(O_start,O_end(3),O_end(4));
O_start(7)= O_start(6);
targetOrbit(7) = 2*pi;

rm0 = orbitalToCar(O_start(1),O_start(2),O_start(3),O_start(4),O_start(5),O_start(6));

ri = orbitalToCar(O_start(1),O_start(2),O_start(3),O_start(4),O_start(5),O_start(6));
rm1 = orbitalToCar(O_start(1),O_start(2),O_start(3),O_start(4),O_start(5),currentOrbit(7));
ra = O_start(1)*(1 + O_start(2));
rm1 = rm1/norm(rm1) * ra;

Orbit1 = orbitfrom2points(ri,rm1,"theta0",0);

% Executing plane and shape change
rpend = orbitalToCar(O_end(1),O_end(2),O_end(3),O_end(4),O_end(5),0);
thetaM3 = 2*pi - acos(dot(rm0,rpend)/(norm(rm0)*norm(rpend)));
rm3 = orbitalToCar(O_end(1),O_end(2),O_end(3),O_end(4),O_end(5),thetaM3);
Orbit2 = orbitfrom2points(rm1,rm3,"theta0",0);
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
orbitDraw(Orbits);
%orbitpropagator(Orbits,15,1000,[1,1,0],"peri");

%% FAST MANOUVER 2.1.3 - Circularize before manouver
% MANOUVER LIST%
% [CHANGE OF SHAPE + PERIAXIS]
% [CHANGE OF PLANE + CHANGE OF SHAPE]
% [CHANGE OF SHAPE + ]
clear
close all
clc

O_start =    [9723.68854055203  0.0865268303717344  0.461525983980154  1.17569878936122  0.982554186272932  1.40668276722613  6.28318530717959];
O_end =      [16720  0.250200000000000  1.11900000000000  0.624500000000000  3.13500000000000  0  3.1];

% In order to save time the first orbit reduces its shape 
[dV,tof,currentOrbit,targetOrbit] = planeChange(O_start,O_end(3),O_end(4));
O_start(7)= O_start(6);
targetOrbit(7) = 2*pi;

rm0 = orbitalToCar(O_start(1),O_start(2),O_start(3),O_start(4),O_start(5),O_start(6));

Orbit1 = O_start;
ri = norm(rm0);
Orbit1(1) = ri;
Orbit1(2) = 0;


rm1 = orbitalToCar(O_start(1),O_start(2),O_start(3),O_start(4),O_start(5),currentOrbit(7));
rm1 = rm1/norm(rm1) * ri;

thetaM1 = acos(dot(rm1,rm0)/(norm(rm0)*norm(rm1)));
Orbit1(7) = O_start(6) + thetaM1;

% Executing plane and shape change
rpend = orbitalToCar(O_end(1),O_end(2),O_end(3),O_end(4),O_end(5),0);
thetaM3 = 2*pi - acos(dot(rm0,rpend)/(norm(rm0)*norm(rpend)));
rm3 = orbitalToCar(O_end(1),O_end(2),O_end(3),O_end(4),O_end(5),thetaM3);
Orbit2 = orbitfrom2points(rm1,rm3,"theta0",0);
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
orbitDraw(Orbits);
%orbitpropagator(Orbits,15,1000,[1,1,0],"peri");

%% FAST MANOUVER 2.1.4 - Circularize at apoapsis before manouver
% MANOUVER LIST%
% [CHANGE OF SHAPE + PERIAXIS]
% [CHANGE OF PLANE + CHANGE OF SHAPE]
% [CHANGE OF SHAPE + ]
clear
close all
clc

O_start =    [9723.68854055203  0.0865268303717344  0.461525983980154  1.17569878936122  0.982554186272932  1.40668276722613  6.28318530717959];
O_end =      [16720  0.250200000000000  1.11900000000000  0.624500000000000  3.13500000000000  0  3.1];

% CIRCULARIZE AT APOAPSIS
O_start(7) = pi;
a0 = O_start(1);
e0 = O_start(0);


% DV CALCULATION
DV1 = simpleDVcalculator(O_start,Orbit1);
DV2 = simpleDVcalculator(Orbit1,Orbit2,1000);
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
orbitDraw(Orbits);
%orbitpropagator(Orbits,15,1000,[1,1,0],"peri");

%% FAST MANOUVER 2.1.3 Cahange periapsis before change
% MANOUVER LIST%
% [CHANGE OF SHAPE + PERIAXIS]
% [CHANGE OF PLANE + CHANGE OF SHAPE]
% [CHANGE OF SHAPE + ]
clear
close all
clc

O_start =    [9723.68854055203  0.0865268303717344  0.461525983980154  1.17569878936122  0.982554186272932  1.40668276722613  6.28318530717959];
O_end =      [16720  0.250200000000000  1.11900000000000  0.624500000000000  3.13500000000000  0  3.1];

% DEFINING FIRST MANOUVER NODE
[dV,tof,currentOrbit,targetOrbit] = planeChange(O_start,O_end(3),O_end(4));
O_start(7)= currentOrbit(7);

rN0 = orbitalToCar(O_start(1),O_start(2),O_start(3),O_start(4),O_start(5),O_start(7));
rA0 = orbitalToCar(O_start(1),O_start(2),O_start(3),O_start(4),O_start(5),pi);
rn0 = orbitalToCar(O_start(1),O_start(2),O_start(3),O_start(4),O_start(5),O_start(6));

Dtheta = acos(dot(rA0,rN0)/(norm(rA0)*norm(rN0)));
omega2 = O_start(5) + Dtheta;
[dV,tof,O_start,Orbit1] = changePeriapsisArg(O_start,omega2);
Orbit1(7) = pi;

% DEFINING SECOND MANOUVER NODE
rpend = orbitalToCar(O_end(1),O_end(2),O_end(3),O_end(4),O_end(5),0);
thetaM3 = 2*pi - acos(dot(rn0,rpend)/(norm(rn0)*norm(rpend)));

% Executing plane and shape change
rm3 = orbitalToCar(O_end(1),O_end(2),O_end(3),O_end(4),O_end(5),0);
rA1 = orbitalToCar(Orbit1(1),Orbit1(2),Orbit1(3),Orbit1(4),Orbit1(5),pi);
Orbit2 = orbitfrom2points(rA1,rm3,"theta0",0);
O_end(6) = 0;

Orbits = [O_start;Orbit1;Orbit2;O_end];

% DV CALCULATION
DV1 = simpleDVcalculator(O_start,Orbit1);
DV2 = simpleDVcalculator(Orbit1,Orbit2);
DV3 = simpleDVcalculator(Orbit2,O_end);

DV_tot = DV1+DV2+DV3;

% TIME CALCULATION
Dt0 = timeOfFlight(O_start);
Dt1 = timeOfFlight(Orbit1);
Dt2 = timeOfFlight(Orbit2);
Dt3 = timeOfFlight(O_end);

Time = Dt0 + Dt1 + Dt2 + Dt3;

fprintf('\nTime optimized manouver 2');
fprintf('\nTotal DV spent is: %f km/s', DV_tot);
fprintf('\nTime in flight: %f s\n', round(Time));

% PLOT
orbitDraw(Orbits);
%orbitpropagator(Orbits,15,1000,[1,1,0],"peri");



%% FAST MANOUVER 2.2 Optimized first transfer
% MANOUVER LIST%
% [CHANGE OF SHAPE + PERIAXIS]
% [CHANGE OF PLANE + CHANGE OF SHAPE]
% [CHANGE OF SHAPE + ]
clear
close all
clc

O_start =    [9723.68854055203  0.0865268303717344  0.461525983980154  1.17569878936122  0.982554186272932  1.40668276722613  6.28318530717959];
O_end =      [16720  0.250200000000000  1.11900000000000  0.624500000000000  3.13500000000000  0  3.1];

% In order to save time the first orbit reduces its shape 
[dV,tof,currentOrbit,targetOrbit] = planeChange(O_start,O_end(3),O_end(4));
O_start(7)= O_start(6);

n = 10;
m = 10;
tol = 0.001;
omega = linspace(0,2*pi,n);
rmanouver1 = linspace(7,14,m) * 1e3;
DV = zeros(m,n);
time = zeros(m,n);
Optimizedorbit = zeros(n,7);

ri = orbitalToCar(O_start(1),O_start(2),O_start(3),O_start(4),O_start(5),O_start(7));
rm1 = orbitalToCar(O_start(1),O_start(2),O_start(3),O_start(4),O_start(5),currentOrbit(7));
rm1Dir = rm1/norm(rm1);
rm3 = orbitalToCar(O_end(1),O_end(2),O_end(3),O_end(4),O_end(5),0);

for j = 1 : m
    executiontime = tic;
    r2 = rm1Dir*rmanouver1(j);
    for i = 1 : n
        Optimizedorbit(i,:) = orbitfrom2points(ri,r2,"omega",omega(i));
    end
    for i = 1 : n
    if Optimizedorbit(i,1) < 0      % Elimino orbite iperboliche
        Optimizedorbit(i,:) = nan;
    end
    r = orbitalToCar(Optimizedorbit(i,1),Optimizedorbit(i,2),Optimizedorbit(i,3),Optimizedorbit(i,4),Optimizedorbit(i,5),Optimizedorbit(i,7));
    if norm(r-r2) > tol
        Optimizedorbit(i,:) = nan;    % Elimino orbite non secanti
    end
    end

    Orbit2 = orbitfrom2points(r2,rm3,"theta0",0);

    for i = 1:n
    if isnan(Optimizedorbit(i,1))
        DV(j,i) = nan;
        time(j,i) = nan;
    else
        DV1 = simpleDVcalculator(O_start,Optimizedorbit(i,:));
        DV2 = simpleDVcalculator(Optimizedorbit(i,:),Orbit2);
        DV(j,i) = DV1 + DV2;
        time(j,i) = timeOfFlight(Optimizedorbit(i,:));
    end

    if time(j,i) > 10000 || DV(j,i) > 10
        DV(j,i) = nan;
        time(j,i) = nan;
    end
    end

    toc(executiontime);
    completion = j/m*100;
    timeleft = (m - j) * executiontime / 60; 
    fprintf('\nLoading... %f percent            ', completion )
    fprintf('Time left: %f min\n', timeleft)
end

%% PLOT
figure Name 'DV optimization'
hold on;
grid on;
ylim([0,2*pi]);
xlim([7,14]*1e3);
zlim([min(min(DV)),max(max(DV))]);
colormap parula;
s = surf(omega,rmanouver1,DV,time);
s.EdgeColor = 'none';

set(gca,'color','k','xcolor','w','ycolor','w','zcolor','w');
set(gcf,'color','k');

ylabel('Radius of manouver node [km]','interpreter','latex','FontSize', 15);
xlabel('omega [rad]','interpreter','latex','FontSize', 15);
zlabel('DV [km/s]','interpreter','latex','FontSize', 15);

axis 'auto xy'

colorbar eastoutside Color 'w' TickLabelInterpreter 'latex'


%% CHOOSING OPTIMAL PARAMETERS

DeltaVopt = min(DV,[],'all');
[r,c] = find(DV == DeltaVopt);
ropt = rmanouver1(r);
omegaopt = omega(c);

r2 = rm1Dir*rmanouver1(j);
ra0 = orbitalToCar(Orbit2(1),Orbit2(2),Orbit2(3),Orbit2(4),Orbit2(5),0);
rpend = orbitalToCar(O_end(1),O_end(2),O_end(3),O_end(4),O_end(5),0);
thetaM3 = acos(dot(ra0,rpend)/(norm(ra0)*norm(rpend)));

Orbit1 = orbitfrom2points(ri,r2,"omega",omegaopt);
Orbit2 = orbitfrom2points(r2,rm3,"theta0",pi);
Orbit2(7) = thetaM3;
O_end(6) = 0;

Orbits = [O_start;Orbit1;Orbit2;O_end];

% DV CALCULATION
DV1 = simpleDVcalculator(O_start,Orbit1,1000);
DV2 = simpleDVcalculator(Orbit1,Orbit2,1000);
DV3 = simpleDVcalculator(Orbit2,O_end,1000);

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
orbitDraw(Orbits);
%orbitpropagator(Orbits,15,1000,[1,0,1],"dynamic",40);

%% FAST MANOUVER 2.3 Leave as it is (Non optimized direct transfer)
% MANOUVER LIST%
% [CHANGE OF SHAPE + PERIAXIS]
% [CHANGE OF PLANE + CHANGE OF SHAPE]
% [CHANGE OF SHAPE + ]
clear
close all
clc

O_start =    [9723.68854055203  0.0865268303717344  0.461525983980154  1.17569878936122  0.982554186272932  1.40668276722613  6.28318530717959];
O_end =      [16720  0.250200000000000  1.11900000000000  0.624500000000000  3.13500000000000  0  3.1];

% In order to save time the first orbit reduces its shape 
[dV,tof,currentOrbit,targetOrbit] = planeChange(O_start,O_end(3),O_end(4));
O_start(7)= currentOrbit(7);
rm1 = orbitalToCar(O_start(1),O_start(2),O_start(3),O_start(4),O_start(5),currentOrbit(7));
Orbit1 = O_start;

% Executing plane and shape change
rm2 = orbitalToCar(O_end(1),O_end(2),O_end(3),O_end(4),O_end(5),0);
Orbit2 = orbitfrom2points(rm1,rm2,"theta0",0);

Orbits = [Orbit1;Orbit2;O_end];

% DV CALCULATION
DV2 = simpleDVcalculator(Orbit1,Orbit2);
DV3 = simpleDVcalculator(Orbit2,O_end);

DV_tot = DV2+DV3;

% TIME CALCULATION
Dt1 = timeOfFlight(Orbit1);
Dt2 = timeOfFlight(Orbit2);
Dt3 = timeOfFlight(O_end);

Time = Dt1 + Dt2 + Dt3;

fprintf('\nTime optimized manouver 2');
fprintf('\nTotal DV spent is: %f km/s', DV_tot);
fprintf('\nTime in flight: %f s\n', round(Time));

% PLOT
orbitDraw(Orbits);
%orbitpropagator(Orbits,15,1000,[1,1,0],"peri");

%% FAST MANOUVER 3 - OPTIMIZATION
% MANOUVER LIST%
% [CHANGE OF SHAPE + PERIAXIS]
% [CHANGE OF PLANE + CHANGE OF SHAPE]
% [CHANGE OF SHAPE + ]
clear
close all
clc

O_start =    [9723.68854055203  0.0865268303717344  0.461525983980154  1.17569878936122  0.982554186272932  1.40668276722613  6.28318530717959];
O_end =      [16720  0.250200000000000  1.11900000000000  0.624500000000000  3.13500000000000  0  3.1];

% In order to save time the first orbit reduces its shape 
[dV,tof,currentOrbit,targetOrbit] = planeChange(O_start,O_end(3),O_end(4));
O_start(7)= O_start(6);
targetOrbit(7) = 2*pi;

rm0 = orbitalToCar(O_start(1),O_start(2),O_start(3),O_start(4),O_start(5),O_start(6));

ri = orbitalToCar(O_start(1),O_start(2),O_start(3),O_start(4),O_start(5),O_start(6));
rm1 = orbitalToCar(O_start(1),O_start(2),O_start(3),O_start(4),O_start(5),currentOrbit(7));
rm1 = rm1/norm(rm1) * 14000;

Orbit1 = orbitfrom2points(ri,rm1,"omega",1.9891);

% Executing plane and shape change
rm2 = orbitalToCar(O_end(1),O_end(2),O_end(3),O_end(4),O_end(5),2*pi);
rpend = orbitalToCar(O_end(1),O_end(2),O_end(3),O_end(4),O_end(5),0);
thetaM3 = 2*pi - acos(dot(rm0,rpend)/(norm(rm0)*norm(rpend)));

n = 200;
m = 200;
tol = 0.001;
omega = linspace(0,2*pi,n);
thetaf = linspace(0,2*pi,m);
DV = zeros(m,n);
time = zeros(m,n);
Optimizedorbit = zeros(n,7);

for j = 1 : m
    executiontime = tic;
    rm3 = orbitalToCar(O_end(1),O_end(2),O_end(3),O_end(4),O_end(5),thetaf(j));
    for i = 1 : n
        Optimizedorbit(i,:) = orbitfrom2points(rm1,rm3,"omega",omega(i));
    end
    for i = 1 : n
    if Optimizedorbit(i,1) < 0      % Elimino orbite iperboliche
        Optimizedorbit(i,:) = nan;
    end
    r = orbitalToCar(Optimizedorbit(i,1),Optimizedorbit(i,2),Optimizedorbit(i,3),Optimizedorbit(i,4),Optimizedorbit(i,5),Optimizedorbit(i,7));
    if norm(r-rm3) > tol
        Optimizedorbit(i,:) = nan;    % Elimino orbite non secanti
    end
    end

    %orbitDraw(Optimizedorbit);

    O_end(6) = thetaf(j);

    parfor i = 1:n
    if isnan(Optimizedorbit(i,1))
        DV(j,i) = nan;
        time(j,i) = nan;
    else
        DV2 = simpleDVcalculator(Orbit1,Optimizedorbit(i,:));
        DV3 = simpleDVcalculator(Optimizedorbit(i,:),O_end);
        DV(j,i) = DV2 + DV3;
        time(j,i) = timeOfFlight(Optimizedorbit(i,:)) + timeOfFlight(O_end);
       
    end

    if time(j,i) > 20000 || DV(j,i) > 10
        DV(j,i) = nan;
        time(j,i) = nan;
    end
    end

    toc(executiontime);
    completion = j/m*100;
    timeleft = (m - j) * executiontime / 60; 
    fprintf('\nLoading... %f percent            ', completion )
    fprintf('Time left: %f min\n', timeleft)
end

%% PLOT
figure Name 'DV optimization'
hold on;
grid on;
xlim([0,2*pi]);
ylim([0,pi/8]);
zlim([min(min(DV)),max(max(DV))-2]);
colormap parula;
s = surf(thetaf,omega,DV,time);
s.EdgeColor = 'none';

set(gca,'color','k','xcolor','w','ycolor','w','zcolor','w');
set(gcf,'color','k');

ylabel('Theta [rad]','interpreter','latex','FontSize', 15);
xlabel('omega [rad]','interpreter','latex','FontSize', 15);
zlabel('DV [km/s]','interpreter','latex','FontSize', 15);

axis 'auto xy'

colorbar eastoutside Color 'w' TickLabelInterpreter 'latex'

%% CHOOSING OPTIMAL PARAMETERS

DeltaVopt = min(DV,[],'all');
[r,c] = find(DV == DeltaVopt);
thetafopt = thetaf(r);
omegaopt = omega(c);
 
O_end(6) = thetafopt;
rm3 = orbitalToCar(O_end(1),O_end(2),O_end(3),O_end(4),O_end(5),thetafopt);
Orbit2 = orbitfrom2points(rm1,rm3,"omega",omegaopt);

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
orbitDraw(Orbits);
%orbitpropagator(Orbits,15,1000,[1,0,1],"dynamic",40);