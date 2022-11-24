function DV = simpleDVcalculator(Orbit1,Orbit2,tol)
% simpleDVcalculator - allow to calculate DV for a manouver between 2
% given Orbits
%
% PROTOTYPE:
% DV = simpleDVcalculator(Orbit1,Orbit2)
%
% DESCRIPTION:
% The function allows to calculate the necessary DV to pass between the 2 orbits from the common
% node
% IMPORTANT: 
% Final position of orbit 1 [Orbit1(7)] must coincide with initial position
% of Orbit 2 [Orbit2(6)];
%
% INPUT:
% Orbit1    [3x1]   Radius vector of first node     [N/D]
% Orbit2    [3x1]   Radius vector of second node    [N/D]
% tol       [1x1]   Tollerance between two manouvers [km] - Default is 0.01
%
% OUTPUT:
% DV        [1x1]   Vector of orbital parameters    [N/D]
% 

[r1,v1] = orbitalToCar(Orbit1(1),Orbit1(2),Orbit1(3),Orbit1(4),Orbit1(5),Orbit1(7));
[r2,v2] = orbitalToCar(Orbit2(1),Orbit2(2),Orbit2(3),Orbit2(4),Orbit2(5),Orbit2(6));

if nargin < 3
    tol = 0.01;
end

if norm(r2-r1) > tol
    error('Error: the two nodes do not coincide');
end

DV = norm(v2-v1);
