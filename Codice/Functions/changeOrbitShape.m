function [Orbit,DV,Manouvers] = changeOrbitShape(Orbit1,Orbit2,impulse)
% Program executes a orbit shape change at the next apo/peri apsis 
%
% \\\\\\\\\\\\\\\\\\\\\\\\\\
%       TO BE EDITED
% \\\\\\\\\\\\\\\\\\\\\\\\\\
%
% PROTOTYPE:
% [Dv1,Dv2,Dv3,thetaf,Dt] = changeOrbitShape(aI,eI,omegaI,aF,eF,omegaF,rb)
%
% DESCRIPTION:
% Allows to condce a tangent manouver at the closest possible point.
% TBA
%
% INPUT:
% 
%
% OUTPUT:
% 

%% NOTE PER CAPA
% SANGENNARO AIUTAMI TU

%% UTILS
if ismac
    load("../Data/utils.mat",'mu');
else
    load("..\Data\utils.mat",'mu');
end

%%  CHECK ON VARIABLES

if impulse > 3
    error("Number of impulses not compatible with the current orbit")
end

%% VARIABLES

%% DECLARATION OF VARIABLES

Orbit = zeros(impulse+1,7);
DV = zeros(impulse+1,1);


a1 = Orbit1(1);
e1 = Orbit1(2);
a2 = Orbit2(1);
e2 = Orbit2(2);

ra1 = a1*(1+e1);
rp1 = a1*(1-e1);
ra2 = a2*(1+e2);
rp2 = a2*(1+e2);

%% 1 IMPUT MANOUVER
% Manouver activates if impulse numeber is 1 or if initial and final orbits
% are tangent in one of the points.


%% 2 Imput manouver - HomannTransfer
if rA1 ~= rA2 % E non ci sono condizioni esterne ri o Dtimposto
    disp('Executing Order 66')
    % Conditions on thetaI:
    if thetaI > 0 && thetaI <= pi
        r1 = rA1;
        deltaTh = pi - thetaI;
        if eI * eF >= 0         % pi to 0
            thetaf = pi;
            r2 = rP2;
            e12 = (r2-r1)/(r1+r2);
            if r1 >= r2
                omega12 = omegaI;
                thetaManouver = pi;
            else
                omega12 = omegaI+pi;
                thetaManouver = 0;
            end
            
            disp('pito0');
        else
            thetaf = 0;         % pi to pi
            r2 = rA2;
            e12 = (r2-r1)/(r1+r2);
            omega12 = omegaI;
            thetaManouver = pi;
            disp('pitopi');
        end
    else
        r1 = rP1;
        deltaTh = 2*pi - thetaI;
        if eI * eF >= 0         % 0 to 0
            r2 = rA2;
            thetaf = 0;
            thetaManouver = pi; 
        else                    % 0 to pi
            r2 = rP2;
            thetaf = pi;
            thetaManouver = pi;
        end       
    end

    % Calculation of DV:
    p12 = r1 + r2;

    Dv1 = sqrt(2*mu*(1/r1 - 1/p12)) - sqrt(2*mu*(1/r1 - 1/pI));
    Dv2 = sqrt(2*mu*(1/r2 - 1/pF)) - sqrt(2*mu*(1/r2 - 1/p12));
    Dv3 = 0;
    Dt = 0;

    % Plot orbit:
    if plotorbit == 1
        kep1 = [aI,eI,iI,OmegaI,omegaI,thetaI];
        kep12 = [p12/2,e12,iI,OmegaI,omega12,thetaManouver];
        
        [X12,Y12,Z12,V12] = plotOrbit(kep12, stepTh, 0, pi);
        if thetaI ~= 0 && thetaI ~= pi
            [X1,Y1,Z1,V1] = plotOrbit(kep1, stepTh, 0, deltaTh);
            X = [X1,X12];
            Y = [Y1,Y12];
            Z = [Z1,Z12];
            V = [V1,V12];
        else
            X = X12;
            Y = Y12;
            Z = Z12;
            V = V12;
        end
    end
end
