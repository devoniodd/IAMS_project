clear
close all
clc

load("../Data/utils.mat",'mu');

e = 0.5;
theta = linspace(0,2*pi,1000);

E_method0 = @(theta) 2 * atan(sqrt((1-e)/(1+e)) * tan(theta/2));
E_method1 = @(theta) acos((e + cos(theta))./(1+e*cos(theta)));
E_method2 = @(thata) asin((sqrt(1-e^2)*sin(theta))./(1+e*cos(theta)));
E_method3 = 2 * atan(sqrt((1-e)/(1+e)) * tan(theta/2));
E_method3(501:end) = E_method3(501:end) + 2*pi;
E_method4 = @(theta) (2 * atan(sqrt((1-e)/(1+e)) * tan(theta/2))).*(theta<pi) + pi.*(theta==pi) +(2*pi + 2 * atan(sqrt((1-e)/(1+e)) * tan(theta/2))).*(theta>pi);

figure
hold on;
plot(theta,E_method0(theta),'r');
plot(theta,E_method1(theta),'b');
plot(theta,E_method2(theta),'g');
%plot(theta,E_method3,'k');
plot(theta,E_method4(theta));

%% Plot di th
clear
close all
clc

e = [0.9,0.8,0.7,0.6,1./(2:9) , 1./(10.^(1 : 10))];
theta = linspace(0,2*pi,1000);
E_method4 = @(theta,e) (2 * atan(sqrt((1-e)/(1+e)) * tan(theta/2))).*(theta<pi) + pi.*(theta==pi) +(2*pi + 2 * atan(sqrt((1-e)/(1+e)) * tan(theta/2))).*(theta>pi);
colormap = winter(22);

figure Name 'Andamento di T in funzione di e'
hold on;
for n = 1 : 10
    plot(theta,E_method4(theta,e(n)),"Color",colormap(n,:));
end