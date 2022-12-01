function orbitpropagator(O,t_step,warp,propOptions,cameraOptions,camerapos)
% orbitpropagator: allows to plot a dynamic orbti and ground trcking
%
% PROTOTYPE:
% orbitpropagator(O,t_step,warpfactor,CameraOptions,Cameraoptions)
%
% DESCRIPTION:
% The function allows to draw a dynamic orbit on a time-dyscrete systeme.
% Possible configurations are in the description.
%
% INPUT:
% O             [nx7]   Orbit matrix                    [N/D]
% t_step        [1x1]   Delta time between wvwry step   [ s ] [Default is 50]
% warp          [1x1]   Warp factor                     [rad] [Default is 1e4]
% propOpptions  [1x3]   Propagator option vector        [see description]
% cameraOptions [1x3]   Camera option vector            [see description]
%
%=========================== PROPAGATOR OPTIONS ===========================
%|                                                                        |
%| Options = [DynamicOrbitpropagation, groundtrack, exportasvideo];       |
%|                                                                        |
%| DynamicOrbitpropagation   -->     Toggle dynamic orbit propagation;    |
%| [0 --> Off]   [1 --> On ]                                              |
%|                                                                        |
%| Groundtrack               -->     Toggle gound tracking;               |
%| [0 --> Off]   [1 --> On ]                                              |
%|                                                                        |
%| Exportasvideo             -->     Allows to export as video;           |
%| [0 --> Off]   [1 -->.AVI]                                              |
%|                                                                        |
%| Default = [1,1,0]                                                      |
%|                                                                        |
%==========================================================================
%
%============================= CAMERA OPTIONS =============================
%| CAMERA SETTING:                                                        |
%| cameraOptions =   "manual"    -->   "Manual" static;                   |
%|               =   "peri"      -->   Static on perifocal system;        |
%|               =   "dynamic"   -->   Rotation around earth;             |
%|               =   "following" -->   Moving with satellite;             |
%|               =   "track"     -->   Always pointing to satellite;      |
%| Default = "manual"                                                     |
%|                                                                        |
%| ADDITIONAL SETTINGS:                                                   |
%| "manual" -->  camerapos  =   [x,y,z] --> Camera position vector [N/D]  |
%|               --> Default = [1,-0.5,0.5]                               |
%| "dynamic"-->  camerapos  =   [angle] --> Camera elevation       [deg]  |
%|               --> Default = 30                                         |
%| "target" -->  camerapos  =   [x,y,z] --> Camera position        [N/D]  |
%|               --> Default = [1,1,1]                                    |
%|                                                                        |
%==========================================================================
%

%================================== DEBUG =================================
showactualwarpfactor = 0;
clouds = 1;
transition = 1;
Width = 1080;
Height = 720;
hideinfo = 1;
hideui = 1;
hidecolormapinvideo = 0;
linewidthplot = 2;
%==========================================================================

%% DEFINING OPTIONS FOR FSOLVE;
options = optimoptions('fsolve','Display','off');

%% CHECK ON IMPUT;

% CHECK ON CAMERA SETTING
if nargin < 6

    if cameraOptions == "manual"
        disp('Warning Camera Position not set, default value will be used')
        camerapos = [1,-0.5,0.5];
    end

    if cameraOptions == "dynamic"
        disp('Warning Camera Angle not set, default value will be used')
        camerapos = 30;
    end

    if cameraOptions == "target"
        disp('Warning Camera Angle not set, default value will be used')
        camerapos = [1,0,0];
    end

    if nargin < 5
        cameraOptions = "manual";
        camerapos = [1,-0.5,0.5];
        if nargin < 4
            propOptions = [1,1,0];
            if nargin < 3
                warp = 10000;
                if nargin < 2
                    t_step = 50;
                end
            end
        end
    end
end

if cameraOptions == "dynamic"
    [sizeimput1,sizeimput2] = size(camerapos);
    if sizeimput1 ~= 1 || sizeimput2 ~= 1 
        error('ERROR: invalid imput for elevation');
    end
elseif nargin == 6
    [sizeimput1,sizeimput2] = size(camerapos);
    if sizeimput1 ~= 1 || sizeimput2 ~= 3
        error('ERROR: imput vector for camerapos has invalid dimensions');
    end
end

[Rows,Col] = size(O);

% Check on imput Matrix O
if Col ~= 7
    error('Matrix O has invalid dimensions');
end




% Check on proOptions;

if length(propOptions) ~= 3
    error('ERROR: propOptions vector has invalid dimensions');
end

propOrbit = propOptions(1);
if propOrbit ~= 0
    if propOrbit ~= 1
        error('ERROR: Invalid value for propOrbit');
    end
end

groundTrack = propOptions(2);
if groundTrack ~= 0
    if groundTrack ~= 1
        error('ERROR: Invalid value for groundTrack');
    end
end

Export = propOptions(3);
capture = 0;
if Export ~= 0
    capture = 1;
    if Export ~= 1
        error('ERROR: Invalid value for Export');
    end
end

%% UTILS IMPORT 
if ismac
    load("../Data/utils.mat",'mu');
else
    load("..\Data\utils.mat",'mu');
end


%==========================================================================
%% ////////////////////// TIME TO COORDDINATES ////////////////////////////
%==========================================================================

%% DELTAT CALCULATION FOR EVERY MANOUVER
% Allows to have a vector of finite elements at the beginning of every manouver
Deltat = zeros(Rows+1,1);
for i = 1 : Rows
    Deltat(i+1) = timeOfFlight(O(i,:),O(i,6),O(i,7));
end

%% DEFINING ABSOLUTE TIMESPAN AND COORDINATE SYSTEM
t = (0 : t_step : sum(Deltat));
CAR = zeros(4,length(t));
index = zeros(1,Rows+1);

%% ORBITALTOCAR FOR EVERY MANOUVER

for i = 1 : Rows
    %% PREALLOCATION OF TEMPORARY MATRIX
    start_time = timeOfFlight(O(i,:),0,O(i,6));
    t_inflight = (start_time : t_step : start_time+Deltat(i+1));
    E_vect = zeros(1,length(t_inflight)+1);
    th_vect = zeros(1,length(t_inflight));
    CAR_t = zeros(4,length(t_inflight));

    %% IMPORTING PARAMETERS
    a = O(i,1);
    e = O(i,2);
    inclination = O(i,3);
    Omega = O(i,4);
    omega = O(i,5);

    %% TIME-->E / E-->THETA DECLARATION OF ANONIMOUS FUNCTIONS

    p = a*(1-e^2);
    timeE = @(E) sqrt(a^3/mu) .* (E - e.*sin(E));       % ELLIPTIC AND CIRCULAR ORBITS
    timeP = @(D) 0.5 * sqrt(p^3 / mu) * (D - D^3 / 3);  % PARABOLICAL ORBITS
    timeH = @(F) sqrt(-a^3 / mu) * (e*sinh(F) - F);     % HYPERBOLIC ORBITS
    
    thE = @(E) (2 * atan(sqrt((1+e)/(1-e)) * tan(E/2))).*(E<pi) + pi.*(E==pi) +(2*pi + 2 * atan(sqrt((1+e)/(1-e)) * tan(E/2))).*(E>pi);
    thP = @(D) 2*atan(D);
    thH = @(F) 2*atan(sqrt((1+e)/(e-1)) * tanh(F/2));

    % CHECK ON ECCENTRICITY
    if  e >= 0 && e < 1
        time = @(E) timeE(E);
        th = @(E) thE(E);
    elseif e == 1
        time = @(D) timeP(D);
        th = @(D) thP(D);
    else 
        time = @(F) timeH(F);
        th = @(F) thH(F);
    end

    %% TIME --> THETA
    for j = 1 : length(t_inflight)                          
        time_it = @(E) time(E) - t_inflight(j);             % DEFINING INTERSECTION QUERY
        E_vect(j+1) = fsolve(time_it,E_vect(j),options);    % TIME              --> ECCENTRIC ANOMALY
        th_vect(j) = th(E_vect(j+1));                       % ECCENTRIC ANOMALY --> THETA
    end
    
    %% ORBITALTOCAR FOR EVERY THETA(t)
    % CAR: Cartesian Coordinates + Speed MATRIX
    % 
    %                  | X(t) |
    % CAR_VECTOR(t) =  | Y(t) |
    %                  | Z(t) |
    %                  | V(t) |
    
    for j = 1 : length(t_inflight)
        [r,v] = orbitalToCar(a,e,inclination,Omega,omega,th_vect(j));
        CAR_t(1,j) = r(1);      % X(t)
        CAR_t(2,j) = r(2);      % Y(t)
        CAR_t(3,j) = r(3);      % Z(t)
        CAR_t(4,j) = norm(v);   % v(t)
    end

    %% IMPORTING COORDINATES IN ABSOLUTE COORDINATES SYSTEM
    cCAR_t = size(CAR_t,2);
    CAR(:,index(i)+1 : index(i)+cCAR_t) = CAR_t;
    index(i+1) = index(i)+cCAR_t;
    cCAR = size(CAR,2);

    
end

%% ADJUSTING TIME
if size(t) ~= cCAR
    t = 0 : t_step : t_step*cCAR;
end

%==========================================================================
%% //////////////////// SETTING UP DYNAMIC PLOTTING ///////////////////////
%==========================================================================

%% SETTING 3D GRADIENT
Vcolorcode = round(CAR(4,:),2)*100;     % Speed is rounded to become a finite number to link to a RGB triplet
Vcolorcode = ceil(Vcolorcode);          % Double check to make every value of speed an integer number;
minV = min(Vcolorcode) - 1;             
maxV = max(Vcolorcode);
Vcolorcode = Vcolorcode - minV;         % Defining speed interval;
cc = jet(maxV-minV);                    % Creates a RGB matrix with specified colormap; 
color = zeros(3,length(CAR(4,:)));      % Preallocation of Color ---> Color(v) = [R(v),G(v),B(v)];
for i = 1:1:length(CAR(4,:))-1          % Color(v) --> Color(v(t))
    color(:,i) = cc(Vcolorcode(i),:)'; 
end

%% DEFINING EARTH
% Creating Earth Shape
nPol = 50;                              % Number of polygons for Earth shape
E_radius = 6378.1363;
E_flattening = 0.0033528131;            % Flattening of Earth ellypsoid;
E_z = E_radius*(1-E_flattening);        % Z size of eath;

% Earth rotation
Erotstep = 360 / (24*60*60) * t_step;   % Defining rotation of Earth for every t_step
zaxys = [0,0,1];
origin = [0,0,0];

%% SETTING ENVIRONMENT FOR ORBIT PROPAGATION
if propOrbit == 1                           % THE FOLLOWING ACTIONS WILL BE EXECUTED ONLY IF ORBITPROP IS ACTIVE

OP = figure(Name = 'Orbit propagation');    % Open figure
OrbitP = axes('Parent',OP);                 % Creates a OrbitProp specific axis sistem
hold on;

% Creating earth shape with textures;
[xE,yE,zE] = ellipsoid(0,0,0,E_radius,E_radius,E_z,nPol);   
Earth = surf(xE,yE,zE);
Esurfce = imread('Earthtexture.jpg');
set(Earth,'Facecolor','texturemap','Edgecolor','none','Cdata',Esurfce);

% setting background to black
set(gca,'color','k','xcolor','w','ycolor','w','zcolor','w');
set(gcf,'color','k');

%Clouds
if clouds == 1
    C_radius = E_radius+20;
    C_z = C_radius*(1-E_flattening);
    [xC,yC,zC] = ellipsoid(0,0,0,C_radius,C_radius,C_z,nPol);
    Clouds = surf(xC,yC,zC);
    Cloudmap = imread('earthcloudmap.jpg');
    Cloudtrans = imread('earthcloudmaptrans.png');
    set(Clouds,'Facecolor','texturemap','Edgecolor','none','Cdata',Cloudmap);
    set(Clouds,'AlphaDataMapping','scaled','AlphaData',Cloudtrans);

    Clouds.FaceColor = 'texturemap'; 
    Clouds.EdgeColor = 'none';
    Clouds.CData = Cloudmap;

    Clouds.FaceAlpha = 'texturemap';
    Clouds.AlphaData = max(Cloudmap,[],3);

    Crotstep = Erotstep*1.2; % Setting rotation of clouds;
end

% FIGURE OPTIONS

grid on;
axis equal;

xlim([min(CAR(1,:)) - 5000,max(CAR(1,:)) + 5000]);
ylim([min(CAR(2,:)) - 5000,max(CAR(2,:)) + 5000]);
zlim([ min([min(CAR(3,:)),E_radius]), max( [max(CAR(3,:)),E_radius] ) ]);

Speed = annotation('textbox', [0, 0.6, 0, 0], 'string', 'Speed:','Color','w');
Tflight = annotation('textbox', [0, 0.5, 0, 0], 'string', 'T:','Color','w');

colormap jet;
caxis([minV/100,maxV/100]);
OPcolorbar = colorbar(OrbitP,'eastoutside',Color='w');
OPcolorbar.Label.String = 'Speed [km/s]';
end

%==========================================================================
%% //////////////////// SETTING UP GROUND TRACKING ////////////////////////
%==========================================================================

TH = zeros(size(CAR(1,:)));     % Preallocation of THETA(t) [Azimuth]   vector
PHI = zeros(size(CAR(1,:)));    % Preallocation of PHI(t)   [Ascension] vector
R = zeros(size(CAR(1,:)));      % Preallocation of R(t)     [Radius] vector
height = zeros(size(CAR(1,:))); % Preallocation of H(t)     [Height]    vector

for i = 1 : length(CAR(1,:))
    [TH(i),PHI(i),R(i)] = cart2sph(CAR(1,i),CAR(2,i),CAR(3,i));    % THETA(X,Y,Z), PHI(X,Y,Z), R(X,Y,Z)
    TH(i) = wrapToPi(TH(i) - deg2rad(Erotstep*(i-1)));             % Setting THETA between -180 & 180
    height(i) = R(i) - E_radius;                                   % H = R(X,Y,Z) - E_radius
end

lat = rad2deg(PHI); % Azimuth to latitude;
lon = rad2deg(TH);  % Elevation to longitude;

if groundTrack == 1
    % Avoids latching lines in graph
    skip = zeros(1,length(lon)-1);
    for i = 1 : length(lon)-1
        if abs(lon(i+1) - lon(i)) > 170
           skip(i) = 1;
        end
    end

    % SETTING GRAPH ENVIRONMENT
    GT = figure(Name = 'Ground Tracking');
    GroundT = axes('Parent',GT);
    hold on;
    axis equal

    xlim([-180,180]);
    ylim([-85,85]);
    xlabel('Longitude','FontSize',15,'Color','k');
    ylabel('Latitude','FontSize',15,'Color','k');
    colormap parula;
    heightmap = parula(length(height));
    
    cdata = imread('earth.png');
    image('CData',flipud(cdata),'XData',[-180,180],'YData',[-90,90])
end

%==========================================================================
%% //////////////////// SETTING UP CAMERA VIEW ////////////////////////////
%==========================================================================
dynamiccamera = 0;
dynamictarget = 0;
if propOrbit
if cameraOptions == "manual"
    maxDistance = max(R);
    camerapos = camerapos/norm(camerapos)*maxDistance;
    OrbitP.CameraPosition = camerapos;
    OrbitP.CameraTarget = [0,0,0];
end

if cameraOptions == "peri"
    dynamiccamera = 1;
    camerapos = zeros(3,length(R));
    maxDistance = max(R);
    for n = 1 : Rows
        it = ceil((index(n+1)- index(n))/3);
        i = O(n,3);
        Omega = O(n,4);
        th = i;
        phi = wrapTo2Pi(Omega-pi/2);
        perif = [cos(phi)*sin(th),sin(phi)*sin(th),cos(th)];
        perif = perif/norm(perif) * maxDistance;
        hold on;
        camerapos(1,index(n)+1 : index(n+1)) = perif(1); 
        camerapos(2,index(n)+1 : index(n+1)) = perif(2);
        camerapos(3,index(n)+1 : index(n+1)) = perif(3);
        if transition == 1 && n < Rows
            nexti = O(n+1,3);
            nextOmega = O(n+1,4);
            nextth = linspace(i,nexti,it);
            nextphi = wrapTo2Pi(linspace(Omega,nextOmega,it)-pi/2);
            transperif = [cos(nextphi).*sin(nextth);sin(nextphi).*sin(nextth);cos(nextth)];
            for j = 1 : it
                transperif(:,j) = transperif(:,j)/norm(transperif(:,j)) * maxDistance;
            end
            camerapos(:,index(n+1)-it+1 : index(n+1)) = transperif;
        end
    end
    OrbitP.CameraTarget = [0,0,0];
end

if cameraOptions == "dynamic"
    maxD = max([ abs(min(CAR(1,:))),abs(min(CAR(2,:))), max(CAR(1,:)), max(CAR(2,:)) ]);
    xlim(OrbitP,[-maxD-5000,maxD+5000]);
    ylim(OrbitP,[-maxD-5000,maxD+5000]);
    dynamiccamera = 1;
    dynamictarget = 0;
    maxDistance = max(R);
    inclination = deg2rad(camerapos);
    camerapos = zeros(3,length(R));
    camerapos(1,:) = cos(wrapTo2Pi(TH'));
    camerapos(2,:) = sin(wrapTo2Pi(TH'));
    camerapos(3,:) = tan(inclination);
    for j = 1 : length(R)
        camerapos(:,j) = camerapos(:,j) / norm(camerapos(:,j)) .* maxDistance;
    end

    OrbitP.CameraTarget = [0,0,0];
end

if cameraOptions == "following"
    dynamiccamera = 1;
    dynamictarget = 0;
    maxDistance = max(R);
    for j = 1 : length(R)
        camerapos(:,j) = CAR(1:3,j) / norm(CAR(1:3,j)) * maxDistance;
    end
    OrbitP.CameraTarget = [0,0,0];
end

if cameraOptions == "target"
    dynamiccamera = 0;
    dynamictarget = 1;
    maxDistance = max(R);
    targetpos = zeros(3,length(R));
    for j = 1 : length(R)
        targetpos(:,j) = CAR(1:3,j);
    end
    OrbitP.CameraPosition = camerapos/norm(camerapos) * maxDistance;
end
end

%==========================================================================
%% ////////////////////////// DYNAMIC PLOT ////////////////////////////////
%==========================================================================
wait = t_step/warp;

%==========================================================================
%/////////////////////////////// VIDEO ////////////////////////////////////
if capture && propOrbit
    OP.Position = [0,0,Width,Height];
    Video = VideoWriter('Orbitprop','Archival');
    %Video.Quality = videoquality;
    Video.FrameRate = ceil(warp/t_step);

    if hideinfo
        Speed.Visible = 'off';
        Tflight.Visible = 'off';
    end
   
    if hideui
        tb = uitoolbar(OP);
        tb.Visible = 'off';
        if hidecolormapinvideo
            OPcolorbar.Visible = 'off';
        end
    end

    wait = 0;

    open(Video);
end
%//////////////////////////////////////////////////////////////////////////
%==========================================================================

tstart = tic;
for i = 1 : length(CAR(1,:)) - 1
    if propOrbit
        animatedline(OrbitP,CAR(1,i:i+1),CAR(2,i:i+1),CAR(3,i:i+1),'Color', color(:,i),'linewidth', linewidthplot);
  
        rotate(Earth,zaxys,Erotstep,origin);
        if clouds == 1
            rotate(Clouds,zaxys,Crotstep,origin);
        end
        set(Speed,'String',"Speed norm:" + CAR(4,i) + "km/s");
        set(Tflight,'String',"Time in flight:" + t(i) + "s");
  
        if dynamiccamera == 1
            OrbitP.CameraPosition = camerapos(:,i);
        end
        if dynamictarget == 1
            OrbitP.CameraTarget = targetpos(:,i);
        end
    end

    if groundTrack && skip(i) == 0
       line(GroundT,lon(i:i+1),lat(i:i+1),'Color',heightmap(i,:));
    end

    if capture
        frame = getframe(OP);     
        writeVideo(Video, frame);
    end

    drawnow limitrate;
    pause(wait);
end
telapsed = toc(tstart);

if showactualwarpfactor == 1
    actualwarp = sum(Deltat)/telapsed;
    fprintf('Actual warp factor is: %d \n',actualwarp);
end

if capture && propOrbit
    close(Video)
    close(OP);
    disp('Video successfully captured and saved');
end

return;
