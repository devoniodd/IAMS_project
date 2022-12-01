delete(gcp('nocreate')); % Close every multi core environment
clear
close all
clc

%% 3 IMPULSE MANOUVER

% This matlab script calcultes an optimized 3 impulse manouver with
% relative DV spent and time of flight.
% The first manouver consist of a circularization at apoapsis, the
% second is a DV optimized direct transfer from the intersection between
% the two orbital planes to a generic node located on the final orbit.
% Please note that the second manouver will change more than 1 orbital
% parameter at a time, therefore for ease of coding (not necesserely for 
% optimization), transfer orbits will be calculated with the function 
% Orbitfrom2point.m, and DV will be obtained vectorially with the function 
% simpleDVcalculator.m 

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

% FINDING MANOUVER NODE
[dV,tof,currentOrbit,targetOrbit] = planeChange(O_start,O_end(3),O_end(4));
manouvernode = orbitalToCar(currentOrbit(1),currentOrbit(2),currentOrbit(3),currentOrbit(4),currentOrbit(5),currentOrbit(6));
manouvernodeDir = manouvernode/norm(manouvernode);

% Obtaining the direction of the line of intersection of the orbital
% planes.

% CIRCULARIZE AT APOAPSIS
O_start(7) = pi;
a0 = O_start(1);
e0 = O_start(2);
rA0 = a0*(1+e0);

Orbit1 = O_start;
Orbit1(1) = rA0;
Orbit1(2) = 0;
Orbit1(6) = pi;

% Executing circularization at apoapsis

% CHANGE PLANE + CHANGE SHAPE
rm1 = orbitalToCar(Orbit1(1),Orbit1(1),Orbit1(3),Orbit1(4),Orbit1(5),Orbit1(6));    % Obtaining first manouver node position
DthetaM2 = acos( dot(rm1,manouvernodeDir) / ( norm(rm1)*norm(manouvernodeDir) ) );   
thetaM2 = pi + DthetaM2;                                                            % Obtaining real anomaly of second manouver node
Orbit1(7) = thetaM2;                                                                % Defining orbit arc
rm2 = orbitalToCar(Orbit1(1),Orbit1(2),Orbit1(3),Orbit1(4),Orbit1(5),Orbit1(7));    % Obtaining second manouver node position

%% TRANSFER OPTIMIZATION
% Transfer orbit is not defined, in the following section a graph with DV
% and transfer time is calculated in function of 2 arbitrary parameters.
% Parameter 1:      thetaf -->  Determinates the position of arrival node on final orbit;
% Parameter 2:      omega  -->  Determinates the shape and size of trnasfer orbit;
% Note: omega has been chosen as parameter only for ease of display on
% following graphs.
% Note 2: abnormal orbits in shape, transfer time or DV required will be
% instantly discarted.

delete(gcp('nocreate')); % Close every multi core environment

%==========================================================================
%///////// FOLLOWING SECTION REQUIRES PARALLEL COMPUTING TOOLBOX //////////
%==========================================================================


%======================== [for loop OPTIONS] ==============================
n = 200;
m = 400;
maxDV = 10;
maxtime = 3*1e4;
tol = 0.001;

parpool(4); % maximum nuber of workers in parfor loop --> parpool(x), x < Core numbers
%==========================================================================

omega = linspace(-pi,pi,n);     % omega discretization
thetaf = linspace(-pi,pi,m);    % thetaf discretization
DV = zeros(m,n);                % Prelocation of memory
time = zeros(m,n);

timeelapsed = tic;

parfor j = 1 : m
    O_endhypotetical = O_end;
    O_endhypotetical(6) = thetaf(j);
    rm3 = orbitalToCar(O_endhypotetical(1),O_endhypotetical(2),O_endhypotetical(3),O_endhypotetical(4),O_endhypotetical(5),O_endhypotetical(6));
    Optimizedorbit = zeros(n,7);

    for i = 1 : n
        Optimizedorbit(i,:) = orbitfrom2points(rm2,rm3,"omega",omega(i));       % Defining transfer orbits
    end

    for i = 1 : n
        if Optimizedorbit(i,1) <= 0                              % Discarting Hyperbolic orbits
            Optimizedorbit(i,:) = nan;
        end

        rs = orbitalToCar(Optimizedorbit(i,1),Optimizedorbit(i,2),Optimizedorbit(i,3),Optimizedorbit(i,4),Optimizedorbit(i,5),Optimizedorbit(i,6));
        re = orbitalToCar(Optimizedorbit(i,1),Optimizedorbit(i,2),Optimizedorbit(i,3),Optimizedorbit(i,4),Optimizedorbit(i,5),Optimizedorbit(i,7));

        if norm(re-rm3) > tol || norm(rs-rm2) > tol
            Optimizedorbit(i,:) = nan;                          % Discarting orbits without intersections in arrival or departure node
            DV(j,i) = nan;
            time(j,i) = nan;
        else

            if ~isnan(Optimizedorbit(i,1))
                DV2 = simpleDVcalculator(Orbit1,Optimizedorbit(i,:));
                DV3 = simpleDVcalculator(Optimizedorbit(i,:),O_endhypotetical);
                DV(j,i) = DV2 + DV3;
                Dt1 = timeOfFlight(Optimizedorbit(i,:));
                Dt2 = timeOfFlight(O_endhypotetical);
                time(j,i) =  Dt1 + Dt2;
            else
                DV(j,i) = nan;
                time(j,i) = nan;
            end

        end

        if time(j,i) > maxtime || DV(j,i) > maxDV 
            DV(j,i) = nan;                                      % Discarting orbits with eccessive time or DV
            time(j,i) = nan;
        end
    end
end

toc(timeelapsed)

delete(gcp('nocreate')); % Close every multi core environment

%% PLOT DV
% PLOTTING THE OBTAINED DV GRAPH
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

xlabel('$\omega_f$ [rad]','interpreter','latex','FontSize', 15);
ylabel('$\theta_t$ [rad]','interpreter','latex','FontSize', 15);
zlabel('$\Delta$V ($\omega_f$,$\theta_t$) [km/s]','interpreter','latex','FontSize', 15);

axis 'auto xy'

colorbar eastoutside Color 'w' TickLabelInterpreter 'latex'

%% PLOT TIME
% PLOTTING THE OBTAINED TIME GRAPH
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

xlabel('$\omega_f$ [rad]','interpreter','latex','FontSize', 15);
ylabel('$\theta_t$ [rad]','interpreter','latex','FontSize', 15);
zlabel('Transfer time t($\omega_f$,$\theta_t$) [s]','interpreter','latex','FontSize', 15);

axis 'auto xy'

colorbar eastoutside Color 'w' TickLabelInterpreter 'latex'

%% CHOOSING OPTIMAL PARAMETERS FROM OPTIMIZATION MATRIX
% It is possible to optimize the transfer manouver exploiting the previou
% DV graph.
% Note: in the following section optimal parameters are obtained from the 
% previous DV matrix.

DVopt = min(DV,[],"all");

[r,c] = find(DV == DVopt);
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

orbitDraw(Orbits);
%orbitpropagator(Orbits,15,1000,[1,0,0],"peri");

%% ANALYTICAL SOLUTION OF MINIMUM DV PROBLEM
% The transfer orbit obtained with the DV matrix is subject to a coarse
% discretization and therefore will not be precise enough while showing
% the most optimized orbit. 
% From the general behaviour of the plotted orbits it is possible to make 
% some observations:
% - The most optimized transfer orbit tends to approximate the
%   inclination of the first orbit.
% - DV decreases as thetaf value gets closer to thetaN, where thetaN is
%   the theta of anti-direction of the line of orbital plane intersection
% - In the previous case DV is optimal when periapsis(omea) --> departure 
%   manouver node
%
% From the previous observations it is possible to conclude that optimal DV
% is obtained with a cahnge of shape to the final intersection node, and
% a final change of inclination, shape and apsis.
% This conclusion is compliant wiht theory as a tangent manouver at periapsi 
% and a plane change at apoapsis resoults in the most efficient configuaration.
% An analytical solution can thereby be defined as follows:

rAf = orbitalToCar(O_end(1),O_end(2),O_end(3),O_end(4),O_end(5),pi);
thetafopt = acos( dot(rAf,rm2)/(norm(rAf)*norm(rm2)) );
O_end(6) = 2*pi - thetafopt;
rp2 = norm(rm2);
ra2 = O_end(1)* (1-O_end(2)^2) / ( 1 + O_end(2)*cos(thetafopt) );
a2 = (ra2+rp2)/2;
e2 = (ra2-rp2)/(2*a2);
Orbit2 = [a2,e2,Orbit1(3),Orbit1(4),Orbit1(5) + Orbit1(7) - Orbit1(6) - pi ,0,pi];

Orbits = [O_start;Orbit1;Orbit2;O_end];

% DV CALCULATION
DV1 = simpleDVcalculator(O_start,Orbit1);
DV2 = simpleDVcalculator(Orbit1,Orbit2);
DV3 = simpleDVcalculator(Orbit2,O_end,120); % This has to be adjusted - is currently broken

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

%% CHOOSING OPTIMAL PARAMETERS - TIME OPTIMIZATION
% From what has been showed previusly the 3 impulse manouver tends to be
% extremely more DV-efficient than the standard manouver, the following is
% an attempt to optimize transfer time, allowing a similar DV expense to
% the standard manouver.

%======================== [for loop OPTIONS] ==============================
DVtreshhold = 6.4;
%==========================================================================

DV_2 = DV + DV1;
time_2 = time;

for j = 1 : m              % Discarting all solutions with a total DV higher than standard manouver
    for i = 1 : n
        if DV_2(j,i) > DVtreshhold
            time_2(j,i) = nan;
            DV_2(j,i) = nan;
        end
    end
end

%% PLOT
figure Name 'time optimization'
hold on;
grid on;
ylim([0,2*pi]);
xlim([0,pi]);
zlim([min(min(time_2)),max(max(time_2))-2]);
colormap winter;
s = surf(omega,thetaf,time_2,DV_2);
s.EdgeColor = 'none';

set(gca,'color','k','xcolor','w','ycolor','w','zcolor','w');
set(gcf,'color','k');

xlabel('$\omega_t$ [rad]','interpreter','latex','FontSize', 15);
ylabel('$\theta_f$ [rad]','interpreter','latex','FontSize', 15);
zlabel('Manoeuver Time t($\omega_t$,$\theta_f$) [s]','interpreter','latex','FontSize', 15);

axis 'auto xy'

colorbar eastoutside Color 'w' TickLabelInterpreter 'latex'

%% CHOOSING OPTIMAL OPTIONS

timeopt = min(time_2,[],'all');
[r,c] = find(time_2 == timeopt);
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

orbitDraw(Orbits);
%orbitpropagator(Orbits,15,1000,[1,0,0],"peri");

