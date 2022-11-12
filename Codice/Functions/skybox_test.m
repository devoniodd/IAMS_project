clear
close all
clc


skybox = imread("milkyway.jpg");
nPolsky = 50;
figure(Name = 'Galaxy');

ax0 = axes;
axis equal;
hold on;

r = 10e10;
set(gca,'Color','k');

xlim([0,r]);
ylim([-r,r]);
zlim([-r,r]);

% Creating skybox
[xS,yS,zS] = sphere(ax0,nPolsky);
xS = xS * r;
yS = yS * r;
zS = zS * r;
Skybox = surf(ax0,xS,yS,zS,'Facecolor','texturemap','Edgecolor','none','Cdata',skybox);

ax0.Clipping = 'on';
ax0.ClippingStyle = '3dbox';
ax0.OuterPosition = [-2,-2,5,5];
ax0.CameraPosition = [0,0,0];