function orbitDraw(O,options)
% orbitDraw: display orbit plot and manouvers
%
% PROTOTYPE:
% orbitDraw(O,options)
%
% DESCRIPTION:
% Plot orbits from matrix O in subcession.
%
% INPUT:
% O             [nx7]   Orbit matrix                    [N/D]
% options       [1xn]   Options vector                  [N/D] 
%
%
%================================ OPTIONS =================================
%|                                                                        |
%| Options = ["OptionName", "Setting"]                                    |
%|                                                                        |
%| SETTINGS:                                                              |
%|                                                                        |
%| "Ocompletion"             -->     Toggle orbit completion              |
%| ["Yes"]   ["No"]                  Default is "Yes"                     |
%|                                                                        |
%| "Export"                  -->     Export picture to specified path     |
%| ["path"]                                                               |
%|                                                                        |
%| "setViasual"              -->     Set viewpoint position               |
%| [ [Azimuth , Elevation] ]                                              |
%|                                                                        |
%| "thStep"                  -->     Set Dtheta for discretization        |
%| [ th_step ]                       Default is 50                        |
%|                                                                        |
%==========================================================================

%% ============================== DEFAULT =================================
showOrbitCompletion = 1;
th_step = 0.01;
%==========================================================================

%% READING IMPUTS
if nargin >= 2
if (-1)^(length(options)) > 0
    for n = 1 : 2 : length(options)
        % SHOW ORBIT COMPLETION

        if options(n) == "Ocompletion"
            if options(n+1) == "Yes" || options(n+1) == "yes"
                showOrbitCompletion = 1;
            elseif options(n+1) == "No" || options(n+1) == "no"
                showOrbitCompletion = 0;
            else
                error('Invalid expression after "Ocomplition" ')
            end
        end
        
        % EXPORT
        if options(n) == "Export"
            path = options(n+1);
        end

        % INITIAL VIEWPOINT
        if options(n) == "setVisual"
            [r,c] = size( options(n+1) );
            if r ~= 1 || c ~= 2
                error('Invalid expression after "setVisual" ')
            end
            viewPoint = options(n+1);
        end

        % THSTEP
        if options(n) == "thStep"
            [r,c] = size( options(n+1) );
            if r ~= 1 || c ~= 1
                error('Invalid expression after "thStep" ')
            end
            th_step = options(n+1);
        end
    end
else
    error('Error, some options were not given');
end
end

%% STARTING FUNCTIONS
[Rows,Col] = size(O);
if Col ~= 7
    error('Matrix O has invalid dimensions');
end

%% UTILS IMPORT 
if ismac
    load("../Data/utils.mat",'mu');
    addpath(genpath("../Functions/Texture/"))
    
else
    load("..\Data\utils.mat",'mu');
    addpath(genpath("..\Functions\Texture\"))
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


