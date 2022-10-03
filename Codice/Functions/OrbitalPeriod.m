function [T] = OrbitalPeriod(r,M)
% The function allows to calculate the orbital period T from the radius.
% Input:
% r - radius of orbit
% M - mass of the planet (If left unspecified erth mass is used)

% To call this function type[ P ] = orbital period ( R )

%% CONSTANTS
G = 6.673e-11;

%% PERIOD
if nargin == 1
    M = 5.98e24;
end

T = sqrt((4*(pi^2)*r^3)/(G*M));

end