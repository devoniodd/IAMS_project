function [r,v] = orbitalToRV(a,e,i,Omega,omega,theta,mu,degrees)

if nargin == 7
    error(sprintf("Please specify if the given angles are:\n0 - Radiants\n1 - Degrees"));
end

%% UTILS
p = a(1-norm(e)^2);

%% POSITION VECTOR
r = [(p*cos(theta))/(), ()/(), ()/()];

%% VELOCITY VECTOR
v


end

