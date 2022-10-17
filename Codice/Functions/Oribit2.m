clear
close all
clc

% Call the Terra_3D Function

kep = [22e3,0.8,15,45,80,180];
timeScaleFactor = 800;



%% Call the plotOrbit function
[X,Y,Z,V] = plotOrbit(kep,deg2rad(0.1),1,2*pi);

%% Plot the 3D satellite orbit
% Setting colormap
ax1 = axes;
Vcolorcode = round(V,2)*100;
Vcolorcode = ceil(Vcolorcode);
maxV = max(Vcolorcode);

cc = jet(maxV);

%Plotting
Terra3d(ax1);
axis equal;

% Define an indefinite plot
h = plot3(nan,nan,nan,'or');

orbit = plot3(X,Y,Z,'LineWidth',1);

%% Define the step animation
step_animation = 5;

% Define the moving point
dt = zeros(1,length(X));
color = zeros(3,length(X));
for i = 1:step_animation:length(X)-1
    ds = norm([X(i),Y(i),Z(i)] - [X(i+step_animation),Y(i+step_animation),Z(i+step_animation)]);
    dt(i) = (ds/V(i))/timeScaleFactor;
    color(:,i) = cc(Vcolorcode(i),:)';
end

for i = 1:step_animation:length(X)-1
    set(h,'XData',X(i),'YData',Y(i),'ZData',Z(i));
    set(orbit,'color', color(:,i));
    drawnow limitrate
    pause(dt(i));
end


