function Orbit = orbitfrom2points(r1,r2,setparameter,var)
% orbitfrom2point - allows to obtain all the orbit kep parameters trought 2
% given radiuses on the orbit.
%
% PROTOTYPE:
% Orbit = orbitfrom2points(r1,theta1,r2,theta2)
%
% DESCRIPTION:
% The function allows to define an orbit from 2 given nodes.
% Angles must be given in radians.
%
% INPUT:
% r1            [3x1]       Radius vector of first node             [km ]
% r2            [3x1]       Radius vector of second node            [km ]
% setparameter  [string]    Parameter to be set
% var           [1x1]       Value of inserted parameter
% 
%
% OUTPUT:
% Orbit     [7x1]   Vector of orbital parameters (orbit arc)        [N/D]
% 

%% UTILS IMPORT
if ismac
    load("../Data/utils.mat",'mu');
else
    load("..\Data\utils.mat",'mu');
end

%% SET theta0
if setparameter == "theta0"
    theta1 = var;
    %% ORBIT ANOMALY
    rcross = cross(r1,r2);
    ismorethanpi = 0;

    if rcross(3) < 0
        rcross = cross(r2,r1);
        ismorethanpi = 1;
        theta2 = theta1 + 2*pi - acos( dot(r1,r2)/(norm(r1)*norm(r2)) );
    else
        theta2 =  theta1 + acos( dot(r1,r2)/(norm(r1)*norm(r2)) );
    end
    rcross = rcross/norm(rcross);

    %% Inclination;
    iDir = [1,0,0];
    jDir = [0,1,0];
    kDir = [0,0,1];

    i = acos(dot(rcross,kDir));

    %% RAAN
    N = cross(kDir,rcross);
    N = N/norm(N);

    if dot(N,jDir) >= 0
        Omega = acos(dot(N,iDir));
    else 
        Omega = 2*pi - acos(dot(N,iDir));
    end

    %% ARGUMENT OF PERIAPSIS
    r1Dir = r1/norm(r1);
    omega = acos(dot(N,r1Dir)) - theta1;

    if ismorethanpi
        omega = 2*pi - omega;
    end

    %% PERIFOCAL PARAMETERS CALCULATION
    r1_norm = norm(r1);
    r2_norm = norm(r2);

    e = abs((r1_norm-r2_norm)/(r1_norm*cos(theta1) - r2_norm*cos(theta2)));
    p = r1_norm*(1+e*cos(theta1));
    a = p/(1-e^2);
    
    Orbit = [a,e,i,Omega,omega,theta1,theta2];
end

%% SET omega
if setparameter == "omega"
    
    omega = var;
    iDir = [1,0,0];
    jDir = [0,1,0];
    kDir = [0,0,1];

    rcross = cross(r1,r2);
    if rcross(3) < 0
        rcross = cross(r2,r1);
        Dtheta = 2*pi - acos( dot(r1,r2)/(norm(r1)*norm(r2)) );
    else
        Dtheta =  acos( dot(r1,r2)/(norm(r1)*norm(r2)) );
    end
    rcross = rcross/norm(rcross);

    %% RAAN
    N = cross(kDir,rcross);
    N = N/norm(N);

    if dot(N,jDir) >= 0
        Omega = acos(dot(N,iDir));
    else 
        Omega = 2*pi - acos(dot(N,iDir));
    end

    %% TEHTA 1 & THETA2
    if r1(3) >= 0
        theta1 = acos(dot(N,r1)/norm(r1)) - omega;
    else
        theta1 = 2*pi - acos(dot(N,r1)/norm(r1)) - omega;
    end
    theta2 = wrapTo2Pi(theta1+Dtheta);

    %% Inclination;

    i = acos(dot(rcross,kDir));

    %% PERIFOCAL PARAMETERS CALCULATION
    r1_norm = norm(r1);
    r2_norm = norm(r2);

    e = abs((r1_norm-r2_norm)/(r1_norm*cos(theta1) - r2_norm*cos(theta2)));
    p = r1_norm*(1+e*cos(theta1));
    a = p/(1-e^2);
    
    Orbit = [a,e,i,Omega,omega,theta1,theta2];
end


return;
