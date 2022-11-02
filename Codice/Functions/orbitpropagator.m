function [lon,TH] = orbitpropagator(O,t_step,warpfactor)

% Imput is a matrix containing 5 Parameters for orbital characterization +
% 2 Parameters with initial and final real anomalies
% Es O = [a,e,i,O,o,thetaI,thetaF];

% How it works
% 1     DT calcultation for every orbit
% 2     Time discretization with dh
% 3     Theta discretization in funcion of time --> Non linear equation
%       (may be useful an EDO)
% 4     Orbit segment plotting --> CartoKep
% 5     Actual orbit propagation
% 6     Ground tracking?

[Rows,Col] = size(O);
if Col ~= 7
    error('Matrix O has invalid dimensions');
end

options = optimoptions('fsolve','Display','off');

%% UTILS IMPORT 
if ismac
    load("../Data/utils.mat",'mu');
else
    load("..\Data\utils.mat",'mu');
end


%% Deltat calculation for each manouvers
% Allows to have a vector of finite elements at the beginning of every manouver
Deltat = zeros(Rows+1,1);
for i = 1 : Rows
    Deltat(i+1) = timeOfFlight(O(i,:),O(i,6),O(i,7));
end

%% Defining an absolute timespan/coordinates
t = (0 : t_step : sum(Deltat));
CAR = [];

%% Cart2Kep for every manouver

for i = 1 : Rows
    t_inflight = (0 : t_step : Deltat(i+1));
    E_vect = zeros(1,length(t_inflight)+1);
    th_vect = zeros(1,length(t_inflight));
    CAR_t = zeros(4,length(t_inflight));
    a = O(i,1);
    e = O(i,2);
    inclination = O(i,3);
    Omega = O(i,4);
    omega = O(i,5);

    %% Defining theta

    % Metodi possibili: 
    % 1 --> Attraverso il calcolo puntuale della anomalia eccentrica Pro:
    % easy, visto a esercitazione, Contro: richiede soluzione di un probelma non
    % lineare ad ogni iterazione --> Non ottimizzato
    
    time = @(E) sqrt(a^3/mu) .* (E - e.*sin(E));
    dtime = @(E) sqrt(a^3/mu) .* (- e.*cos(E));

    for j = 1 : length(t_inflight)
        time_it = @(E) time(E) - t_inflight(j);
        E_vect(j+1) = fsolve(time_it,E_vect(j),options);
        % E_vect(j+1) = newton(E_vect(j),1000,1e-6,time_it,dtime);
        th_vect(j) = 2*atan( sqrt((1+e)/(1-e)) * tan(E_vect(j+1)/2) );
    end

    %% Kep2Car for every theta
    for j = 1 : length(t_inflight)
        [r,v] = orbitalToCar(a,e,inclination,Omega,omega,th_vect(j));
        CAR_t(1,j) = r(1);
        CAR_t(2,j) = r(2);
        CAR_t(3,j) = r(3);
        CAR_t(4,j) = norm(v);
    end

    CAR = [CAR,CAR_t];
end

% Setting 3D gradient
Vcolorcode = round(CAR(4,:),2)*100; % Approssima le velocita' durante l'orbita e le trasforma in un numero finito
Vcolorcode = ceil(Vcolorcode); % Messa per sicurezza, rende tutti gli elementi del vettore interi
maxV = max(Vcolorcode);
cc = jet(maxV); % Crea una matrice con un grediente di colori RGB
color = zeros(3,length(CAR(4,:))); % Il vettore conterra' i colori associati alla velocita' nel plot dell'orbita
for i = 1:1:length(CAR(4,:))-1 % Associa il colore alla velocita nell'orbita
    color(:,i) = cc(Vcolorcode(i),:)'; 
end

%% Setting Environment


% Creating Earth Shape
nPol = 50;
E_radius = 6378.1363;
E_flattening = 0.0033528131;
E_z = E_radius*(1-E_flattening);

figure(Name = 'Orbit propagation');
hold on;

[xE,yE,zE] = ellipsoid(0,0,0,E_radius,E_radius,E_z,nPol);
Earth = surf(xE,yE,zE);
Esurfce = imread('Earthtexture.jpg');
set(Earth,'Facecolor','texturemap','Edgecolor','none','Cdata',Esurfce);

% setting background to black
set(gca,'color','k','xcolor','w','ycolor','w','zcolor','w');
set(gcf,'color','k');

% Earth rotation

Erotstep = 360 / (24*60*60) * t_step;
zaxys = [0,0,1];
origin = [0,0,0];

%Clouds
C_radius = E_radius+20;
C_z = C_radius*(1-E_flattening);
[xC,yC,zC] = ellipsoid(0,0,0,C_radius,C_radius,C_z,nPol);
Clouds = surf(xC,yC,zC);
Cloudmap = imread('earthcloudmap.jpg');
Cloudtrans = imread('earthcloudmaptrans.png');
set(Clouds,'Facecolor','texturemap','Edgecolor','none','Cdata',Cloudmap);
set(Clouds,'AlphaDataMapping','scaled',...
    'AlphaData',Cloudtrans);

Clouds.FaceColor = 'texturemap'; 
Clouds.EdgeColor = 'none';
Clouds.CData = Cloudmap;

Clouds.FaceAlpha = 'texturemap';
Clouds.AlphaData = max(Cloudmap,[],3);


Crotstep = Erotstep*1.2;

% defining boundaries

grid on;
axis equal;

xlim([min(CAR(1,:)),max(CAR(1,:))]);
ylim([min(CAR(2,:)),max(CAR(2,:))]);
%zlim([min(CAR(3,:)),max(CAR(3,:))]);

Speed = annotation('textbox', [0, 0.6, 0, 0], 'string', 'Speed:','Color','w');
Tflight = annotation('textbox', [0, 0.5, 0, 0], 'string', 'T:','Color','w');




%% Propagazione dell'orbita;
tstart = tic;
for i = 1 : length(CAR(1,:)) - 1
  animatedline(CAR(1,i:i+1),CAR(2,i:i+1),CAR(3,i:i+1),'Color', color(:,i));
  
  rotate(Earth,zaxys,Erotstep,origin);
  rotate(Clouds,zaxys,Crotstep,origin);
  set(Speed,'String',"Speed norm:" + CAR(4,i) + "km/s");
  set(Tflight,'String',"Time in flight:" + t(i) + "s");
  drawnow
  pause(t_step/warpfactor);
end
telapsed = toc(tstart);
actualwarp = sum(Deltat)/telapsed;
fprintf('Actual warp factor is: %d \n',actualwarp);

%% Ground tracking; - da sistemare.
TH = zeros(size(CAR(1,:)));
PHI = zeros(size(CAR(1,:)));
height = zeros(size(CAR(1,:)));

for i = 1 : length(CAR(1,:))
    [TH(i),PHI(i),R] = cart2sph(CAR(1,i),CAR(2,i),CAR(3,i));
    TH(i) = TH(i) - deg2rad(Erotstep*(i-1));
    height(i) = R - E_radius;
end
lat = rad2deg(PHI); % Cambio da elevazione a azimuth
lon = rad2deg(TH);

% Faccio apparire tutto sullo stesssa faccia:
for i = 1 : length(CAR(1,:))
    if lon(i) > 180
        lon(i) = lon(i) - 360*(floor(lon(i)/360));    
    end
    if lon(i) < -180
        lon(i) = lon(i) + 360*(floor(lon(i)/360));
    end
end

figure(Name="Ground tracking");
nlon = [-180,180];
nlat =[-85,85];
geolimits(nlat,nlon);
hold on;
geoscatter(lat,lon,1,height,'.');
colormap turbo;
geobasemap colorterrain;
