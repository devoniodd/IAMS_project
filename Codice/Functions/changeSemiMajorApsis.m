function [dV,Orbit1,Orbit2] = changeSemiMajorApsis(Orbit1,afinal)
% Simple 1 impulse manouver function that changes orbit dimension.

%% UTILS
if ismac
    load("../Data/utils.mat",'mu');
else
    load("..\Data\utils.mat",'mu');
end

%% DECLARATION OF VARIABLES
a1 = Orbit1(1);
e1 = Orbit1(2);
omegaOld = Orbit1(5);
theta0 = Orbit1(6);

%% DV CALCULATION
rP1 = a1*(1-e1);
rA1 = a1*(1+e1);

if theta0 > pi && theta0 < 2*pi
    r1 = rP1;
    thetaM = 2*pi;
else 
    r1 = rA1;
    thetaM = pi;
end

r2 = 2*afinal - r1;
efinal = abs((r1-r2)/(2*afinal));

if r2 < 6400
    error('New radius is too close to earth surface');
end

if thetaM == 2*pi
    if r2 > r1
        omegaNew = omegaOld;
        thetafinal = 0;
    else
        omegaNew = wrapToPi(omegaOld + pi);
        thetafinal = pi;
    end
else
    if r2 > r1
        omegaNew = wrapToPi(omegaOld + pi);
        thetafinal = 0;
    else
        omegaNew = omegaOld;
        thetafinal = pi;
    end
end
p1 = a1*(1-e1^2);
p2 = afinal*(1-efinal^2);
dV = sqrt(mu/p2)*(1-efinal*cos(thetaM)) - sqrt(mu/p1)*(1-e1*cos(thetafinal));

Orbit1(7) = thetaM;
Orbit2 = [afinal,efinal,Orbit1(3),Orbit1(4),omegaNew,thetafinal,wrapTo2Pi(thetafinal+pi)];