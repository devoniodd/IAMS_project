%% ES.3b
clear all;
clc;
%% Inputs
a=7000;
mu=398600;
e=0.7108;
t0=0;
P=(2.000001)*pi*sqrt(a^3/mu);
f0=0;
e=[0,0.2,0.4,0.4,0.6,0.8,0.95];
t=linspace(t0,2*P,200);
%% Computing of true anomaly for different eccentricity values
for j=1:length(e)
for i=1:length(t)
y(i)=trueanomaly(a,e(j),mu,t0,f0,t(i));
end
f(:,j)= y;
end
f=f';
%% Plots
figure(1)
plot(t,f(1,:),'r')
hold on
plot(t,f(2,:),'y')
plot(t,f(3,:),'b')
plot(t,f(4,:),'g')
plot(t,f(5,:),'k')
plot(t,f(6,:),'g')
plot(t,f(7,:),'c')
grid on;
xlabel('time [s]')
ylabel('true anomaly [rad]')
pbaspect([1.25 1 1])
legend('e=0', 'e=0.2', 'e=0.4', 'e=0.6', 'e=0.8', 'e=0.95')

figure(2)
[X,Y]=meshgrid(e,t);
Z=f';
surf(Y,X,Z)
xlabel('time [s]')
ylabel('eccentricity [s]')
zlabel('true anomaly [rad]')

