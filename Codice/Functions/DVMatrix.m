delete(gcp('nocreate')); % Close every multi core environment
clear
close all
clc

%% DV OPTIMIZED DIRECT


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

%% BRUTEFORCE APPROACH
delete(gcp('nocreate')); % Close any multi core environment
%======================== [for loop OPTIONS] ==============================
n = 100;
m = 100;
o = 200;
tol = 0.001;
hmax = 200;

%parpool(6); % maximum nuber of workers in parfor loop --> parpool(x), x < Core numbers
%==========================================================================

Rearth = 6378;
Rmax = Rearth + hmax;
thetaD = linspace(0,2*pi,n);
thetaA = linspace(0,2*pi,m);
omega = linspace(0,2*pi,o);
Points = zeros(n,m,o);
DV = zeros(n,m,o);
omegapos = zeros(n,m,o);

telapsed = tic;

parfor i = 1 : n % CYCLE ON THETAI
    thetai = thetaD(i);
    [ri,vi] = orbitalToCar(aS,eS,iS,OmegaS,omegaS,thetai);
    for j = 1 : m % CYCLE ON THETAF
        thetaf = thetaA(j);
        [rf,vf] = orbitalToCar(O_end(1),O_end(2),O_end(3),O_end(4),O_end(5),thetaf);

         for k =  1 : o % CYCLE ON ORBIT SHAPE
            orbitT = orbitfrom2points(ri,rf,"omega",omega(k));
            
            [r1,v1] = orbitalToCar(orbitT(1),orbitT(2),orbitT(3),orbitT(4),orbitT(5),orbitT(6));
            [r2,v2] = orbitalToCar(orbitT(1),orbitT(2),orbitT(3),orbitT(4),orbitT(5),orbitT(7));

            if orbitT(1) > 0 && orbitT(2) > 0
                if orbitT(1)*(1 - orbitT(2)) < Rmax
                    DV(i,j,k) = nan;
                    omegapos(i,j,k) = nan;
                end
                if norm(r1-ri) < tol && norm(r2-rf) < tol
                    DV1 = norm(vi-v1);
                    DV2 = norm(v2-vf);
                    DV(i,j,k) = DV1 + DV2;
                    omegapos(i,j,k) = omega(k);
                else
                    DV(i,j,k) = nan;
                    omegapos(i,j,k) = nan;
                end
            else
                DV(i,j,k) = nan;
                omegapos(i,j,k) = nan;
            end
        end

    end

end

toc(telapsed)

delete(gcp('nocreate')); % Close any multi core environment

%% MIN DV

[v,loc] = min(DV(:));
[x,y,z] = ind2sub(size(DV),loc);

thetaDopt = thetaD(x);
thetaAopt = thetaA(y);
omegaopt = omega(z);

O_start(7) = thetaDopt;
O_end(6) = thetaDopt;

[rD,vD] = orbitalToCar(O_start(1),O_start(2),O_start(3),O_start(4),O_start(5),O_start(7));
[rA,vA] = orbitalToCar(O_end(1),O_end(2),O_end(3),O_end(4),O_end(5),O_end(6));

Orbitt = orbitfrom2points(rD,rA,"omega",omegaopt);

Orbits = [O_start;Orbitt;O_end];

% DV CALCULATION
DV1 = simpleDVcalculator(O_start,Orbitt);
DV2 = simpleDVcalculator(Orbitt,O_end);

DV_tot = DV1+DV2;

% TIME CALCULATION
Dt0 = timeOfFlight(O_start);
Dt1 = timeOfFlight(Orbitt);
Dt2 = timeOfFlight(O_end);

Time = Dt0 + Dt1 + Dt2;

fprintf('\nTime optimized manouver 2');
fprintf('\nTotal DV spent is: %f km/s', DV_tot);
fprintf('\nTime in flight: %f s\n', round(Time));

orbitDraw(Orbits);
%orbitpropagator(Orbits,15,1000,[1,0,0],"peri");

%% PLOT DV

% PLOTTING THE OBTAINED DV GRAPH
figure Name 'DV optimization'
hold on;
grid on;
colormap parula;
set(gca,'color','k','xcolor','w','ycolor','w','zcolor','w');
set(gcf,'color','k');
axis equal
colorbar eastoutside Color 'w' TickLabelInterpreter 'latex'

xlabel('$\theta_i$ [rad]','interpreter','latex','FontSize', 15);
ylabel('$\theta_f$ [rad]','interpreter','latex','FontSize', 15);
zlabel('$\omega_t$ [rad]','interpreter','latex','FontSize', 15);

for k = 1 : o
    surf(thetaD,thetaA,omegapos(:,:,k),DV(:,:,k),"EdgeColor","none");
end