function [h,i,N,Omega,e,omega,theta] = rvToOrbital(r,v,mu,degrees)
%% VERSORS
iDir = [1; 0; 0];
jDir = [0; 1; 0];
kDir = [0; 0; 1];

%% VECTOR NORMALIZATION
rNorm = norm(r);
vNorm = norm(v);
vRadial = r*v/rNorm;

%% ANGULAR MOMENTUM
h = cross(r,v);
hNorm = sqrt(h*h);

%% INCLINATION
i = acos(hNorm*kDir/hNorm);

%% NODAL PLANE
N = cross(kDir,h);
NNorm = norm(N);

%% RIGHT ASCENSION OF THE ASCENDING NODE
if NNorm*jDir >= 0
    Omega = acos(NNorm*iDir./NNorm);
else 
    Omega = 360 - acos(NNorm*iDir./NNorm);
end

%% ECCENTRICITY
e = 1/mu * [cross(v,h) - mu*r./rNorm];
eNorm = norm(e);

%% ARGUMENT OF PERIGEE
omega = acos()

%% RADINATS TO DEGREES
if nargin == 3 || degrees = 1 

end
-

end