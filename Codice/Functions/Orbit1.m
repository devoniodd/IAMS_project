clear
close all
clc

% Call the Terra_3D Function

kep = [22e3,0.8,15,45,30,180];
timeScaleFactor = 400;

%% Call the plotOrbit function
[X,Y,Z,V] = plotOrbit(kep,deg2rad(0.1),1,4*pi);
%% Plot the 3D satellite orbit
% Setting colormap
min = round(min(V));
max = round(max(V));

%Plotting
ax1 = axes;
Terra3d(ax1);
axis equal;

ax2 = axes;

% Define an indefinite plot
h = plot3(ax2,nan,nan,nan,'or');

patch(ax2, X,Y,Z,V,'Facecolor','none','Edgecolor','interp');

% Link two axes together
hLink = linkprop([ax1,ax2],{'XLim','YLim','ZLim','CameraUpVector','CameraPosition','CameraTarget'});

% Hide the top axes
ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];

%% Define the step animation
step_animation = 1;
% Define the moving point
dt = zeros(1,length(X));
for i = 1:step_animation:length(X)-1
    ds = norm([X(i),Y(i),Z(i)] - [X(i+step_animation),Y(i+step_animation),Z(i+step_animation)]);
    dt(i) = (ds/V(i))/timeScaleFactor;
end

for i = 1:step_animation:length(X)-1
    set(h,'XData',X(i),'YData',Y(i),'ZData',Z(i));
    drawnow limitrate
    pause(dt(i));
end
