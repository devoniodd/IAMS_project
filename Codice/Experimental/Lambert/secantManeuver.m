function [dV,tof,phiV,orb1,orb2] =  secantManeuver(orb1,orb2,domega,opt)
% secantManeuver - Change of the periapsis anomaly and shape through a secant maneuver
%
% PROTOTYPE:
% [dv,tof,ph1V,orb1,orb2] = secantManeuver(orb1,orb2,domega,opt)
%
% DESCRIPTION:
% Change of the Periapsis anomaly and shape by giving the function the new
% orbital parameters and the omega difference through a secant maneuver.
% Omega Difference is considered positive if counterclockwise.
% Input 0 in optimization to minimize the time required, input 1 to minimize
% the delta v required for the maneuver. 
% (Note that in some cases the result might be the same)
%
% INPUT:
% initial orbit     [1x7]           Initial orbital parameters                          [N/D]
% final orbit       [1x7]           Final orbital parameters                            [N/D]
% delta omega       [1x1]           Difference between anomalies                        [rad]
% optimization      [1x1]           Choice to prioritize either time or delta V         [-]
%
% OUTPUT:
% delta V           [1x1]           Difference between velocities                       [km/s]
% time of flight    [1x1]           Time to tranfer from orbit 1 to orbit 2             [s]
% ph1 V             [1x1]           Flight path angle                                   [rad]
% initial orbit     [1x7]           Initial orbital parameters                          [N/D]
% final orbit       [1x7]           Final orbital parameters                            [N/D]

%% UTILS
if ismac
    load("../../Data/utils.mat",'mu');
else
    load("..\..\Data\utils.mat",'mu');
end

%% VARIABLES CHECK
if domega > 2*pi || domega < 0
    domega = wrapTo2Pi(domega);
end

if opt ~= 1 && opt ~= 0
    error("The optimization parameter must be either 0 or 1")
end

if orb1(3) ~= orb2(3)
    error("The two orbits aren't coplanar and the maneuver can't be applied")
end

if orb1(4) ~= orb2(4)
    error("The two orbits aren't coplanar and the maneuver can't be applied")
end

%% INPUT EXTRACTION

% Eccentricity
e1 = orb1(2);
e2 = orb2(2);

% Semilatus rectum
p1 = orb1(1)*(1-e1^2);
p2 = orb2(1)*(1-e2^2);

% Angular momentum
h1 = sqrt(p1 * mu);
h2 = sqrt(p2 * mu);

%% INTERSECTION POINTS

% Problem definition
a = e1*h2^2 - e2*h1^2 * cos(domega);
b = -e2*h1^2 * sin(domega);
c = h1^2 - h2^2;
phi = atan(b/a);

% Intersection point anomalies
theta1 = wrapTo2Pi(phi + acos(c/a * cos(phi)));
theta2 = wrapTo2Pi(phi - acos(c/a * cos(phi)));

%% MANOUVER POINT

if opt == 1
    % Maneuver point anomaly
    if theta1>theta2
        thetaMP = theta2;
    else 
        thetaMP = theta1;
    end
    
    orb1(7) = thetaMP;
    orb2(5) = wrapTo2Pi(orb1(1,5) + domega);
    orb2(6) = wrapTo2Pi(thetaMP - domega);
    orb2(7) = nan;

    % Radius
    r = p1 * 1/(1+e1*cos(thetaMP));
    
    % Velocities
    Vt1 = h1/r;
    Vr1 = mu/h1 * e1 * sin(thetaMP);
    V1 = sqrt(Vt1^2 + Vr1^2);
    
    Vt2 = h2/r;
    Vr2 = mu/h2 * e2 * sin(orb2(6));
    V2 = sqrt(Vt2^2 + Vr2^2);
    
    % Flight path angle
    gamma1 = atan(Vr1/Vt1);
    gamma2 = atan(Vr2/Vt2);
end

if opt == 0
     % Maneuver point anomaly
    if orb1(6) > theta1 && orb1(6) < theta2
        thetaMP = theta2;
    else 
        thetaMP = theta1;
    end

    orb1(1,7) = thetaMP;
    orb2(1,5) = wrapTo2Pi(orb1(5) + domega);
    orb2(1,6) = wrapTo2Pi(thetaMP - domega);
    orb2(1,7) = nan;
    
    % Radius
    r = p1 * 1/(1+e1*cos(thetaMP));
    
    % Velocities
    Vt1 = h1/r;
    Vr1 = mu/h1 * e1 * sin(thetaMP);
    V1 = sqrt(Vt1^2 + Vr1^2);
    
    Vt2 = h2/r;
    Vr2 = mu/h2 * e2 * sin(orb2(6));
    V2 = sqrt(Vt2^2 + Vr2^2);
    
    % Flight path angle
    gamma1 = atan(Vr1/Vt1);
    gamma2 = atan(Vr2/Vt2);
end
%% OUTPUTS

% Velocities
dV = sqrt(V1^2 + V2^2 - 2*V1*V2*cos(gamma2-gamma1));

% Time of Flight
tof = timeOfFlight(orb1,orb1(1,6),orb1(1,7));

% Flight path Angle
phiV = wrapTo2Pi(atan((Vr2-Vr1)/(Vt2-Vt1)));







