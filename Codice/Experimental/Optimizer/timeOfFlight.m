function [t] = timeOfFlight(orbit,theta1,theta2)
% timeOfFlight - Time required to go from theta1 to theta2
%
% PROTOTYPE:
% [t] = timeOfFlight(orbit,theta1,theta2)
%
% DESCRIPTION:
% The function calculates the time neeeded to tread between two points from
% an orbit of any shape.
% Angles must be given in radians.
%
% INPUT:
% orbit     [1x7]   Orbital parameters vector       [N/D]
% theta1    [1x1]   Starting point of your route    [rad]
% theta2    [1x1]   Ending point of your route      [rad]
%
% OUTPUT:
% t         [1x1]   Time of flight                  [s]
% 
%% UTILS IMPORT
if ismac
    load("../../Data/utils.mat",'mu');
else
    load("..\..\Data\utils.mat",'mu');
end

%% INPUT EXTRACTION
a = orbit(1,1);
e = orbit(1,2);

%% CIRCULAR ORBIT
if  e == 0

    T = 2*pi*sqrt((a^3)/(mu));
    dTheta = wrapTo2Pi(theta2-theta1);
    t = T * dTheta/(2*pi);
    return;

end

%% ELLIPTICAL ORBIT
if e > 0 && e < 1
   
    T = 2*pi/sqrt(mu) * a^(3/2);

    E1 = 2 * atan(sqrt((1-e)/(1+e)) * tan(theta1/2));
    E2 = 2 * atan(sqrt((1-e)/(1+e)) * tan(theta2/2));

    t1 = sqrt(a^3 / mu) * (E1 - e*sin(E1));
    t2 = sqrt(a^3 / mu) * (E2 - e*sin(E2));

    if t2 >= t1
        t = t2 - t1;
    else
        t = t2 - t1 + T;
    end
    return;

end

%% PARABOLICAL ORBIT
if e == 1

    rp = a*(1-e);
    p = 2*rp;

    D1 = tan(theta1/2);
    D2 = tan(theta2/2);

    t1 = 0.5 * sqrt(p^3 / mu) * (D1 - D1^3 / 3);
    t2 = 0.5 * sqrt(p^3 / mu) * (D2 - D2^3 / 3);

    t = abs(t2 - t1);
    return;

end

%% HYPERBOLICAL ORBIT
if e > 1
    
    F1 = 2 * atanh(sqrt((1+e)/(e-1)) * tan(theta1/2));
    F2 = 2 * atanh(sqrt((1+e)/(e-1)) * tan(theta2/2));

    t1 = sqrt(-a^3 / mu) * (e*sinh(F1) - F1);
    t2 = sqrt(-a^3 / mu) * (e*sinh(F2) - F2);

    t = abs(t2 - t1);
    return;

end

end
