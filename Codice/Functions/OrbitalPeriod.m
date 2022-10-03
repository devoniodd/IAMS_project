function [T] = orbitalPeriod(a,M)
% The function allows to calculate the orbital period T from the radius.
% INPUT:
% a - ellipse major apsis
% M - mass of the planet (If left unspecified erth mass is used)
%
% OUTPUT:
% T - orbital period

% To call this function type[ P ] = orbital period ( R )

%% CONSTANTS
G = 6.673e-11;

%% PERIOD
if nargin == 1
    M = 5.98e24;
end

T = sqrt((4*(pi^2)*a^3)/(G*M));

end