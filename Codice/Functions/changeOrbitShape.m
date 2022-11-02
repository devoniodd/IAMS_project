function [Dv1,Dv2,Dv3,thetaf,Dt,X,Y,Z,V] = changeOrbitShape(aI,eI,iI,omegaI,OmegaI,thetaI,aF,eF,Plot,rb,Dtin,Time,Fuel)
% Homann's bitang manouver
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
% FINO A PLOT POSSO USARE NARGIN
% DA RB DEVO USARE EXIST
% SANGENNARO AIUTAMI TU

%% UTILS
if ismac
    load("../Data/utils.mat",'mu');
else
    load("..\Data\utils.mat",'mu');
end

%% Check on variables

if Plot ~= 1 && Plot ~= 0
    error("Please insert a valid value for plot:")
end

if  Plot == 1
    plotorbit = 1;
end

%% Calculating a for transfer orbit

pI = 2*aI; % Major apsis for I orbit
pF = 2*aF; % Major apsis for II orbit

rP1 = aI * (1 - eI); % Periapsis of initial Orbit
rA1 = aI * (1 + eI); % Apoapsis of initial Orbit
rP2 = aF * (1 - eF); % Periapsis of final Orbit
rA2 = aF * (1 + eF); % Apoapsis of final Orbit

%% ThetaF:

%% Default values
stepTh = 0.01;

%% 1 Imput manouver
if rA1 == rA2
    disp("Executing 1 input manouver - simple change of orbit");
    % Conditions on thetaI
    if thetaI > 0 && thetaI <= pi
        r = rA1;
        deltaTh = pi - thetaI;
        if eI * eF >= 0
            thetaf = 0;
        else
            thetaf = pi;
        end
    else
        r = rP1;
        deltaTh = 2*pi - thetaI;
        if eI * eF >= 0
            thetaf = pi;
        else
            thetaf = 0;
        end       
    end

    % DeltaV calculation:
    Dv1 = sqrt(2*mu*(1/r - 1/pF)) - sqrt(2*mu*(1/r - 1/pI)); % Check plz
    
    % Plot orbit:
    if plotorbit == 1
        kep = [pF,eF,iI,OmegaI,omegaf,pi];
        [X,Y,Z,V] = plotOrbit(kep, stepTh, 0, deltaTh, thetaI);
    end
end

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
