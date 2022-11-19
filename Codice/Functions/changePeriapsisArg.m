function [dV,tof,orb1,orb2] = changePeriapsisArg(orb1,domega)
% changePeriapsisArg - Change of the Periapsis anomaly  (Secant Maneuver)
%
% PROTOTYPE:
% [dv,tof,orb1,orb2] = changePeriapsisArg(orb1,domega)
%
% DESCRIPTION:
% Change of the Periapsis anomaly by giving the function the delta omega required.
% Omega Difference is considered positive if counterclockwise. 
%
% INPUT:
% initial orbit     [1x7]           Initial orbital parameters                          [N/D]
% delta omega       [1x1]           Difference between anomalies                        [rad]
%
% OUTPUT:
% delta V           [1x1]           Difference between velocities                       [km/s]
% time of flight    [1x1]           Time to tranfer from orbit 1 to orbit 2             [s]
% initial orbit     [1x7]           Initial orbital parameters                          [N/D]
% final orbit       [1x7]           Final orbital parameters                            [N/D]

%% UTILS
if ismac
    load("../Data/utils.mat",'mu');
else
    load("..\Data\utils.mat",'mu');
end
%% VARIABLES CHECK
if domega > 2*pi || domega < 0
    domega = wrapTo2Pi(domega);
end

%% FINAL ORBIT

% Vector creation
orb2 = zeros(1,7);

% Unchanged Parameters
orb2(1,1) = orb1(1,1);
orb2(1,2) = orb1(1,2);
orb2(1,3) = orb1(1,3);
orb2(1,4) = orb1(1,4);

% Final omega
orb2(1,5) = wrapTo2Pi(orb1(1,5) + domega);

% Final theta
if orb1(1,6) < domega/2 || orb1(1,6) > domega/2 + pi
    orb1(1,7) = domega/2;
else
    orb1(1,7) = do/2 + pi;
end

% Initial and Final theta new orbit
if orb1(1,6) < domega/2 || orb1(1,6) > domega/2 + pi
    orb2(1,6) = 2*pi - domega/2;
else
    orb2(1,6) = pi - domega/2;
end

orb2(1,7) = nan;

%% VELOCITY DIFFERENCE
dV = 2*sqrt(mu/(orb1(1,1)*(1-orb1(1,2)^2)))*orb1(1,2)*sin(domega/2);

%% TIME OF FLIGHT
tof = timeOfFlight(orb1,orb1(1,6),orb1(1,7));


end