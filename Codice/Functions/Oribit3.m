clear
close all
clc

% Call the Terra_3D Function

kep = [10e3,0.2,15,45,80,180];

% O = [10e3,0.2,deg2rad(15),pi/4,deg2rad(80),0,pi];

timeScaleFactor = 800;

% Define the step animation
step_animation = 5;

% Call the plotOrbit function
[Dv1,Dv2,Dv3,thetaf,Dt,X,Y,Z,V] = changeOrbitShape(12e3,0.4,deg2rad(15),deg2rad(45),deg2rad(80),pi/4,22e3,0.2,1);
%[X,Y,Z,V] = plotOrbit(kep, 0.01, 1, pi);

% Setting 3D gradient
Vcolorcode = round(V,2)*100; % Approssima le velocita' durante l'orbita e le trasforma in un numero finito
Vcolorcode = ceil(Vcolorcode); % Messa per sicurezza, rende tutti gli elementi del vettore interi
maxV = max(Vcolorcode);
cc = jet(maxV); % Crea una matrice con un grediente di colori RGB
color = zeros(3,length(X)); % Il vettore conterra' i colori associati alla velocita' nel plot dell'orbita
for i = 1:1:length(X)-1 % Associa il colore alla velocita nell'orbita
    color(:,i) = cc(Vcolorcode(i),:)'; 
end

% Loading Earth
ax1 = axes;
Terra3d;

% Define an indefinite plot (SATELLITE)
h = plot3(nan,nan,nan,'or');

% Plotting orbit
for i = 1 : length(X) - 1
  line('XData', X(i:i+1),'YData',Y(i:i+1),'ZData',Z(i:i+1),'Color', color(:,i));
end


% Define the moving point
% Il primo ciclo for calcola le posizioni e i delta tempo;
dt = zeros(1,length(X));

% Utilizzo i nodi di CGL per permettere uno spaziamento piu' equo degli
% intervalli nel plot:
iter = 1:step_animation:(length(X) - step_animation);
it = length(iter);

for i = 1 : step_animation : (length(X)-step_animation)
    ds = norm([X(i),Y(i),Z(i)] - [X(i+step_animation),Y(i+step_animation),Z(i+step_animation)]);
    dt(i) = (ds/V(i))/timeScaleFactor;
end

% Il secondo ciclo for permette di mostrare il plot dinamico;
for i = 1:step_animation:length(X) - step_animation
    set(h,'XData',X(i),'YData',Y(i),'ZData',Z(i),'color',color(:,i)');
    line('XData', X(i:i+1),'YData',Y(i:i+1),'ZData',Z(i:i+1),'Color', color(:,i));
    drawnow limitrate
    pause(dt(i));
end


