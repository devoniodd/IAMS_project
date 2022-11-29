function orbitDraw(O,options,path)
% orbitDraw: display orbit plot and manouvers
%
% PROTOTYPE:
% orbitpropagator( ... )
%
% DESCRIPTION:
% The function allows to draw the orbit.
%
% INPUT:
% O             [nx7]   Orbit matrix                    [N/D]
% path 
%

%================================== DEBUG =================================
showOrbitCompletion = 1;
displaymanouvernodes = 1;
th_step = 0.01;
%==========================================================================

[Rows,Col] = size(O);
if Col ~= 7
    error('Matrix O has invalid dimensions');
end

%% UTILS IMPORT 
if ismac
    load("../Data/utils.mat",'mu');
else
    load("..\Data\utils.mat",'mu');
end

%% Defining discretization in theta:

DTheta = sum(O(:,7) - O(:,6));
Th = (0 : th_step : DTheta);

fig = figure('Name','Orbit Plot');
hold on;
axis equal;
grid on;

%% Setting Environment
% Creating Earth Shape
nPol = 50;
E_radius = 6378.1363;
E_flattening = 0.0033528131;
E_z = E_radius*(1-E_flattening);
% Earth rotation
[xE,yE,zE] = ellipsoid(0,0,0,E_radius,E_radius,E_z,nPol);
Earth = surf(xE,yE,zE);
Esurfce = imread('Earthtexture.jpg');
set(Earth,'Facecolor','texturemap','Edgecolor','none','Cdata',Esurfce);

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

% Setting color 
color = jet(Rows);

% setting background to black
set(gca,'color','k','xcolor','w','ycolor','w','zcolor','w');
set(gcf,'color','k');


%% Orbit plot:
for i = 1 : Rows
    a = O(i,1);
    e = O(i,2);
    inc = O(i,3);
    Omega = O(i,4);
    omega = O(i,5);
    thi = O(i,6);
    thf = O(i,7);
    % Discretization in th
    if thf >= thi
        Orbit_th = (thi : th_step : thf);
        Orbit_path_th = (thf : th_step : thi + 2*pi);
    else
        Orbit_th = (thi : th_step : thf + 2*pi);
        Orbit_path_th = (thf : th_step : thi);
    end

    Orbit = zeros(3,length(Orbit_th));
    Orbit_path = zeros(3,length(Orbit_path_th));
    
    for j = 1 : length(Orbit_th)
        [R,V] = orbitalToCar(a,e,inc,Omega,omega,Orbit_th(j));
        Orbit(:,j) = R;
    end

    manouver = Orbit(:,end);
    h = plot3(nan,nan,nan,'.',MarkerSize=22,Color=color(i,:));

    for j = 1 : length(Orbit_path_th)
        [R,V] = orbitalToCar(a,e,inc,Omega,omega,Orbit_path_th(j));
        Orbit_path(:,j) = R;
    end
    plot3(Orbit(1,:),Orbit(2,:),Orbit(3,:),'LineWidth',1.5,LineStyle="-",Color=color(i,:));
    
    set(h,'XData',manouver(1),'YData',manouver(2),'ZData',manouver(3));

    if e < 1 && showOrbitCompletion ~= 0
    plot3(Orbit_path(1,:),Orbit_path(2,:),Orbit_path(3,:),LineStyle='--',Color=color(i,:));
    end
end

if exist("viewPoint","var")
    view(viewPoint);
end

if exist("path","var")
    set(gcf,'Color',[0 0 0]);
    exportgraphics(gcf,strcat(path,'.png'),"ContentType","image","BackgroundColor","black")
    saveas(gcf,path,'fig');
end


