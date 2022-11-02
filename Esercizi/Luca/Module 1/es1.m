gamma= c/(2*m);
w0=sqrt(k/m);
w=sqrt(w0^2-gamma^2);
%% Analitycal solution of damped armonic oscillator
pos= @(t) exp(-gamma*t).*(y0(1).*cos(w*t)+((y0(2)+gamma.*y0(1))/w).*sin(w*t));
vel= @(t) exp(-gamma*t).*(y0(2).*cos(w*t)-((w0^2.*y0(1)+gamma.*y0(2))/w).*sin(w*t));
%% Numerical solution of damped harmonic oscillator
tspan=linspace(0,10,100);
options= odeset('RelTol',1e-13,'AbsTol',1e-14);
[T,Y]=ode113( @(t,y) ode_harmonic_oscill(t,y,w0,gamma),tspan,y0,options);
%% Plots
figure(1)
plot(T,Y(:,1),'-o');
hold on;
fplot(pos,[0,10],'-')
xlabel('Time [T]');
ylabel('Position [L]');
title('Position');
xlim([0 4])

figure(2)
plot(T,Y(:,2),'-o');
hold on;
fplot(vel,[0,10],'-')
xlabel('Time [T]');
ylabel('Velocity [L/T]');
title('Velocity');
xlim([0 4]);