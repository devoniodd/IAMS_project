function [r,v] = OrbitToRv(a,e,i,Omega,omega,theta,degrees)
% Description - TBA
% Default is in radians

%% VALUE CHECK
if nargin < 6
    error("Please insert a valid amount of variables");
end

if nargin == 7 && degrees ~= 1
    if degrees ~= 0
        error("Please select a valid option:    1 - Input in degrees    0 - Input in radians");
    end
    disp("Input is in radians");
else
    disp("Input is in degrees");
end

%% DEG TO RAD
if degrees == 1
    Omega = deg2rad(Omega);
    omega = deg2rad(omega);
    theta = deg2rad(theta);
    i = deg2rad(i);
end

%% RADIUS
