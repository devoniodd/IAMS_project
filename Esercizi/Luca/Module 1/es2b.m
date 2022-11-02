%% Exercise 2b: Integrate numerically a Keplerian perturbed ordit
%Input data: r0 (initial position vector)
%            v0 (initial velocity vector)
%            mu (planetary constant of m2)
%            R  (radius of m2)
%            n number of orbits
%            J2 second zonal harmonic
%Function used: ode_orbit2BP.m
% Contributor: Vincenzo Paolella
% Date: 25/09/2020
%% Elaboration of initial conditions and computing of orbit period
r0=[-7128.137,0,0];
v0= [0,-9.7810,0];
mu=398600;
R=6.3781e+03;
n=10;
J2=0.00108263;
y0=[r0 v0]; %vector 1x6 of initial conditions%
h0=[cross(r0,v0)]; %vector 1x3 of initial angular momentum%
nh0=norm(h0); %magnitude of h0%
nr0=norm(r0);%magnitude of r0%
e0=[(1/mu).*cross(v0,h0)-r0./nr0]; %initial eccentricity vector 1x3%
ne0=norm(e0); %magnitude of e0%
a=(nh0^2/mu)*(1/(1-ne0^2)); %semi-major axis of the elliptical orbit%
P=2*pi*sqrt(a^3/mu); %orbit period%
%% Numerical integtration of 2BP equation of motion
tf=n*P;
steps=500*n;
tspan=linspace(0,tf,steps); %steps of integration%
options=odeset('RelTol',1e-13,'AbsTol',1e-14);
[T,F]=ode113( @(t,f) ode_orbit2BP(t,f,mu),tspan,y0,options);
%F is a 300x6 matrix, its m-row represents the state of vector f (see the odefunction) at the
%m-th time step
for i=1:length(tspan)
    r(i,:)=[F(i,1:3)]; %matrix of position vector, m-th row represents the m-th state%
    v(i,:)=[F(i,4:6)];%matrix of velcoity vector, m-th row represents the m-th state%
    nr(i)=norm(r(i,:)); %vector 1x300 of magnitude of r at m-th time step%
    h(i,:)=cross(r(i,:),v(i,:));%matrix of angular momentum,m-th row represents the m-th state%
    nh(i)=norm(h(i,:));%vector 1x300 of h magnitude at m-th time step"
    e(i,:)=[(1/mu).*cross(v(i,:),h(i,:))-r(i,:)./nr(i)]; %matrix 300x3 of eccentricity vector at m-th state%
    ne(i)=norm(e(i,:)); %vector 1x300 of e magnitude at m-th time step"
    Eps(i)=(norm(v(i,:))^2/2)-mu/nr(i); %specific mechanical energy%
    Doteh(i)=dot(e(i,:),h(i,:)); %dot product e-h)%
    cosf(i)=dot(r(i,:),e(i,:))/(nr(i)*ne(i)); %cosine of true anomaly f%
    vr(i)=dot(v(i,:),r(i,:))/nr(i);
    vtheta(i)=(mu/nh(i))*(1+ne(i)*cosf(i));
end
%% Numerical integtration of perturbed 2BP equation of motion
tspan=linspace(0,tf,steps); 
options=odeset('RelTol',1e-13,'AbsTol',1e-14);
[L,M]=ode113( @(t,f) ode_orbit2BP_perturbed(t,f,mu,J2,R),tspan,y0,options);
for i=1:length(tspan)
    pr(i,:)=[M(i,1:3)]; 
    pv(i,:)=[M(i,4:6)];
    pnr(i)=norm(pr(i,:));
    ph(i,:)=cross(pr(i,:),pv(i,:));
    pnh(i)=norm(ph(i,:));
    pe(i,:)=[(1/mu).*cross(pv(i,:),ph(i,:))-pr(i,:)./pnr(i)];  
    pne(i)=norm(pe(i,:)); 
    pEps(i)=(norm(pv(i,:))^2/2)-mu/pnr(i); 
    pDoteh(i)=dot(pe(i,:),ph(i,:)); 
    pcosf(i)=dot(pr(i,:),pe(i,:))/(pnr(i)*pne(i));
    pvr(i)=dot(pv(i,:),pr(i,:))/pnr(i);
    pvtheta(i)=(mu/pnh(i))*(1+pne(i)*pcosf(i));
end
%% Plots
%Plot of the orbit%
figure(1)
plot3(F(:,1),F(:,2),F(:,3),'o')
grid on;
hold on
plot3(M(:,1),M(:,2),M(:,3),'o')
[x,y,z]=sphere(100);
surf(x,y,z);
x=x*R;
y=y*R;
z=z*R;
surf(x,y,z);
xlabel('r_x [Km]');
ylabel('r_y [Km]');
zlabel('r_z [Km]')
axis equal
title('Quasi-circular medium Earth orbit');

%% Plot of the angular momentum
figure(2)
plot(tspan,nh,'-y');
hold on;
plot(tspan,pnh,'-r');
xlabel('t [s]');
ylabel('||h|| [km^2/s]');
title('Angular momentum') 
grid on;

%% Plot of eccentricity vector
figure(3)
plot(tspan,ne,'-y');
hold on;
plot(tspan,pne,'-r')
xlabel('t [s]')
ylabel('||e||')
title('Eccentricity vector')
grid on;

%% Plot of the dot product e-h
figure(4)
plot(tspan,Doteh,'-y')
hold on;
plot(tspan,pDoteh,'-r')
grid on;
xlabel('t [s]')
title('e-h dot product')

%Plot of specific mechanical energy%
figure(5)
plot(tspan,Eps,'-y')
hold on;
plot(tspan,pEps,'-r')
ylabel('\epsilon [km^2/s^2]')
xlabel('t [s]')
grid on;
title('Specific Energy')

%% Plot  of transversal and radial velocity
figure(6)
plot(tspan,vr,'-y')
hold on;
plot(tspan,vtheta,'-y')
plot(tspan,pvtheta,'*r')
plot(tspan,pvr,'xr')
grid on;
xlabel('t [s]')
ylabel('v_r,v_\theta [Km/s]')
title('Radial and transveral velocity')