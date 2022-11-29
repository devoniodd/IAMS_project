%% FAST MANOUVER 2.1.4 - Circularize at apoapsis before manouver

clear
close all
clc

O_start =    [9723.68854055203  0.0865268303717344  0.461525983980154  1.17569878936122  0.982554186272932  1.40668276722613  6.28318530717959];
O_end =      [16720  0.250200000000000  1.11900000000000  0.624500000000000  3.13500000000000  0  3.1];

% FINDING MANOUVER NODE
[dV,tof,currentOrbit,targetOrbit] = planeChange(O_start,O_end(3),O_end(4));
manouvernode = orbitalToCar(currentOrbit(1),currentOrbit(2),currentOrbit(3),currentOrbit(4),currentOrbit(5),currentOrbit(6));
manouvernodeDir = manouvernode/norm(manouvernode);

% CIRCULARIZE AT APOAPSIS
O_start(7) = pi;
a0 = O_start(1);
e0 = O_start(2);
rA0 = a0*(1+e0);

Orbit1 = O_start;
Orbit1(1) = rA0;
Orbit1(2) = 0;
Orbit1(6) = pi;

% CHANGE PLANE + CHANGE SHAPE
rm1 = orbitalToCar(Orbit1(1),Orbit1(1),Orbit1(3),Orbit1(4),Orbit1(5),Orbit1(6));
DthetaM2 = acos( dot(rm1,manouvernodeDir) / ( norm(rm1)*norm(manouvernodeDir) ) );
thetaM2 = pi + DthetaM2;
Orbit1(7) = thetaM2;

% TRANSFER
rm2 = orbitalToCar(Orbit1(1),Orbit1(2),Orbit1(3),Orbit1(4),Orbit1(5),Orbit1(7));

n = 200;
m = 400;
tol = 0.001;
omega = linspace(-pi,pi,n);
thetaf = linspace(-pi,pi,m);
DV = zeros(m,n);
time = zeros(m,n);
Optimizedorbit = zeros(n,7);

for j = 1 : m
    executiontime = tic;
    rm3 = orbitalToCar(O_end(1),O_end(2),O_end(3),O_end(4),O_end(5),thetaf(j));
    for i = 1 : n
        Optimizedorbitit = orbitfrom2points(rm2,rm3,"omega",omega(i));
    end

    for i = 1 : n

        if Optimizedorbit(i,1) < 0      % Elimino orbite iperboliche
            Optimizedorbit(i,:) = nan;
        end

        rs = orbitalToCar(Optimizedorbit(i,1),Optimizedorbit(i,2),Optimizedorbit(i,3),Optimizedorbit(i,4),Optimizedorbit(i,5),Optimizedorbit(i,6));
        re = orbitalToCar(Optimizedorbit(i,1),Optimizedorbit(i,2),Optimizedorbit(i,3),Optimizedorbit(i,4),Optimizedorbit(i,5),Optimizedorbit(i,7));

        if norm(re-rm3) > tol
            Optimizedorbit(i,:) = nan;    % Elimino orbite non secanti
        end
        if norm(rs-rm2) > tol
            Optimizedorbit(i,:) = nan;    % Elimino orbite non secanti
        end
    end

    O_end(6) = thetaf(j);

    for i = 1:n
    if isnan(Optimizedorbit(i,1))
        DV(j,i) = nan;
        time(j,i) = nan;
    else
        DV2 = simpleDVcalculator(Orbit1,Optimizedorbit(i,:));
        DV3 = simpleDVcalculator(Optimizedorbit(i,:),O_end);
        DV(j,i) = DV2 + DV3;
        time(j,i) = timeOfFlight(Optimizedorbit(i,:)) + timeOfFlight(O_end);
       
    end

    if time(j,i) > 30000 || DV(j,i) > 10
        DV(j,i) = nan;
        time(j,i) = nan;
    end
    end

    toc(executiontime);
    completion = j/m*100;
    itleft = m-j;
    timeleft = itleft * executiontime / 60; 
    fprintf('\nLoading... %f percent            ', completion )
    fprintf('Time left: %f min\n', timeleft)
end

%% PLOT DV
figure Name 'DV optimization'
hold on;
grid on;
ylim([0,2*pi]);
xlim([0,pi]);
zlim([min(min(DV)),max(max(DV))-2]);
colormap parula;
s = surf(omega,thetaf,DV,time);
s.EdgeColor = 'none';

set(gca,'color','k','xcolor','w','ycolor','w','zcolor','w');
set(gcf,'color','k');

ylabel('Theta [rad]','interpreter','latex','FontSize', 15);
xlabel('omega [rad]','interpreter','latex','FontSize', 15);
zlabel('DV [km/s]','interpreter','latex','FontSize', 15);

axis 'auto xy'

colorbar eastoutside Color 'w' TickLabelInterpreter 'latex'

%% PLOT TIME
figure Name 'Time optimization'
hold on;
grid on;
ylim([0,2*pi]);
xlim([0,pi]);
zlim([min(min(time)),max(max(time))-2]);
colormap hot;
s = surf(omega,thetaf,time,DV);
s.EdgeColor = 'none';

set(gca,'color','k','xcolor','w','ycolor','w','zcolor','w');
set(gcf,'color','k');

ylabel('Theta [rad]','interpreter','latex','FontSize', 15);
xlabel('omega [rad]','interpreter','latex','FontSize', 15);
zlabel('Time [s]','interpreter','latex','FontSize', 15);

axis 'auto xy'

colorbar eastoutside Color 'w' TickLabelInterpreter 'latex'

%% CHOOSING OPTIMAL PARAMETERS - DELTA V OPTIMIZED

DeltaVopt = min(DV,[],'all');
[r,c] = find(DV == DeltaVopt);
thetafopt = thetaf(r);
omegaopt = omega(c);
 
O_end(6) = thetafopt;
rm3 = orbitalToCar(O_end(1),O_end(2),O_end(3),O_end(4),O_end(5),thetafopt);
Orbit2 = orbitfrom2points(rm2,rm3,"omega",omegaopt);

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
%orbitpropagator(Orbits,15,1000,[1,0,1],"dynamic",40);

%% CHOOSING OPTIMAL PARAMETERS - DELTA V OPTIMIZED

for j = 1 : m
    for i = 1 : n
        if DV(j,i) > 6
            time(j,i) = nan;
            DV(j,i) = nan;
        end
    end
end

%% PLOT
figure Name 'time optimization'
hold on;
grid on;
ylim([0,2*pi]);
xlim([0,pi]);
zlim([min(min(time)),max(max(time))-2]);
colormap winter;
s = surf(omega,thetaf,time,DV);
s.EdgeColor = 'none';

set(gca,'color','k','xcolor','w','ycolor','w','zcolor','w');
set(gcf,'color','k');

ylabel('Theta [rad]','interpreter','latex','FontSize', 15);
xlabel('omega [rad]','interpreter','latex','FontSize', 15);
zlabel('DV [km/s]','interpreter','latex','FontSize', 15);

axis 'auto xy'

colorbar eastoutside Color 'w' TickLabelInterpreter 'latex'

timeopt = min(time,[],'all');

[r,c] = find(time == timeopt);
thetafopt = thetaf(r);
omegaopt = omega(c);
 
O_end(6) = thetafopt;
rm3 = orbitalToCar(O_end(1),O_end(2),O_end(3),O_end(4),O_end(5),thetafopt);
Orbit2 = orbitfrom2points(rm2,rm3,"omega",omegaopt);

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

%% PLOT
orbitDraw(Orbits);
%orbitpropagator(Orbits,15,1000,[1,0,0],"peri");

