function [t] = timeOfFlight(a,e,theta1,theta2,mass,degrees)
% The function allows to calculate the orbital period T from the radius.

% INPUT:
% a         [1x1][km]           Ellipse major apsis                                                 
% e         [1x1][/]            Eccentricity
% theta1    [1x1][rad][deg]     Initial true anomaly
% theta2    [1x1][rad][deg]     Final true anomaly
% mass      [1x1][kg]           Mass of the planet (If left unspecified or 0 earth mass is used)    
% degrees   [1x1][/]            1 for input angle in degrees
%
% OUTPUT:
% t - time of flight

% To call this function type[ P ] = orbital period ( R )

%% UTILS IMPORT

if ismac
    load("../Data/utils.mat",'G','mEarth');
else
    load("..\Data\utils.mat",'G','mEarth');
end

%% DEG TO RAD
if exist("degrees","var") && degrees == 1
    theta1 = deg2rad(theta1);
    theta2 = deg2rad(theta2);
end

%% PRELIMINARY

if ~exist("mass","var") || mass == 0
    mass = mEarth;
end

mu = G*mass;

%% CIRCULAR ORBIT

if ~exist("e","var") || e == 0

    T = 2*pi*sqrt((a^3)/(mu));
    t = T;

    if ~exist("theta1","var")
        return;
    end

    dTheta = abs(theta1-theta2);
    t = T * dTheta/(2*pi);
    return;
end

%% ELLIPTICAL ORBIT

T = 2*pi/sqrt(mu) * a^(3/2);

E1 = 2 * atan(sqrt((1-e)/(1+e)) * tan(theta1/2));
E2 = 2 * atan(sqrt((1-e)/(1+e)) * tan(theta2/2));

t1 = sqrt(a^3 / mu) * (E1 - e*sin(E1));
t2 = sqrt(a^3 / mu) * (E2 - e*sin(E2));

if theta2 >= theta1
    t = t2 - t1;
else
    t = t2 - t1 + T;
end

if t <= 0
    t = T + t;
end