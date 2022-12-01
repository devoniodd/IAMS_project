clear all;
close all;
clc;

%% PLOT 

plotData = 1;

%% INPUTS

addpath(genpath("../Data/"))
load("utils.mat",'mu');
load("PForbs.mat");

%% INITIAL AND FINAL ORBIT DEFINITION

[r1, V1] = orbitalToCar(orb1(1),orb1(2),orb1(3),orb1(4),orb1(5),orb1(6));

[r2, V2] = orbitalToCar(orb2(1),orb2(2),orb2(3),orb2(4),orb2(5),orb2(6));

%% LAMBERT PARAMETERS

% Starting time of flight
t = 4000;

% Exit parameters
hMin = 200;

% Utils
r0 = 6378;
orbitT = zeros(1,7);
it = 1;


while it == 1 || a*(1-e) < r0+hMin

% Transfer orbit velocity
[Vt1,Vt2] = lambertFunc(r1,r2,t,1);

% Find needed orbital parameters
[a, e] = carToOrbital(r1,Vt1);

% Update
it = it+1;
t = t+10;

end

dVlambert = norm(Vt1-V1);
dVfinal = norm(Vt2-V2);
dVtot = dVlambert + dVfinal;

dVtotList = [dVtot];
ts = [t];

it = 1;

while it == 1 || dVtotList(it) < dVtotList(it-1)
    
    % Update
    t = t+10;
    it = it+1;

    % Transfer orbit velocity
    [Vt1,Vt2] = lambertFunc(r1,r2,t,1);
    
    % Velocities
    dVlambert = norm(Vt1-V1);
    dVfinal = norm(Vt2-V2);
    dVtot = dVlambert + dVfinal;
    dVtotList = [dVtotList dVtot];
    ts = [ts t];

end

%% CHOOSE ONE ORBIT TO DISPLAY

% Set time of flight
vLim = 7;
idx = find(dVtotList<=vLim,1,'first');
tChosen = ts(idx);

% Find corresponding transfer orbit
[Vt1,Vt2] = lambertFunc(r1,r2,t,1);

% Export orbit
[orbitT(1), orbitT(2), orbitT(3), orbitT(4), orbitT(5), orbitT(6)] = carToOrbital(r1,Vt1);
[~, ~, ~, ~, ~, orbitT(7)] = carToOrbital(r2,Vt2);
orb1(7) = orb1(6);
orb2(7) = orb2(6);

% Orbits matrix
orbits = [orb1; orbitT; orb2];

%% PLOT 
if plotData

    f1 = figure();
    chart = plot(ts,dVtotList,'lineWidth',2);
    grid on;
    xlabel("\textbf{time of flight}",'Interpreter','latex','FontSize',15)
    ylabel("\textbf{total $\Delta V$}",'Interpreter','latex','FontSize',15)
    dt = datatip(chart,tChosen,dVtotList(idx),'SnapToDataVertex','on');
    hold off;

    orbitDraw(orbits,[52 15],'../../Images/Lambert/lambertOrbits');
    
    % Save
    saveas(f1,"../../Images/Lambert/lambert",'png');

end

